"""
Discrete-event simulation of DetTreeChain (proof-of-work model) under network delay.

Design rules implemented (main.tex Sec. III-D, III-E, III-G):
  * N lanes; slot q(b) = N*q(parent) + lane (lanes 1..N), genesis slot 0.
  * Each node keeps its own accepted set.  PQP = accepted blocks with slot >= q(c),
    ordered by slot; c = first entry, n = second entry.
  * Window rule (Alg. 1): parent in {c, n}, slot empty, lane link == current lane tip.
  * Retirement: if parent == n -> c retires (n becomes c); elif c has N children -> c retires.
  * Miner of lane a: targets slot N*q(c)+a if empty, otherwise N*q(n)+a.  Memoryless retarget.
  * First-seen wins inside a node (Alg. 1 is applied in arrival order).
  * Fork choice ("longest tree"): every block carries the miner's tree at mining time
    (as N lane tips + its current parent c).  A receiving node adopts that tree if it has
    strictly more blocks than its own (ties: keep own), downloading any blocks it lacks
    (headers-first style), then re-applies its own recently received / displaced blocks
    under Alg. 1.  A tree is represented by its N lane tips: because lane links connect
    every accepted block of a lane, the tree = union of the N lane chains.

Network: full mesh, every (block, node) delay i.i.d.: fixed d or Exp(mean d).
Mining: M miners, M/N per lane, equal hash power, total block rate lam (Poisson).

Baselines are the same code with N=1 (single chain; Alg. 1 reduces to "extend the tip").
"""
import heapq
import math
import random
from bisect import bisect_left, bisect_right, insort


class Sim:
    def __init__(self, N, M, lam, delay, delay_kind, n_blocks, seed,
                 pending_cap=64, sample_every=4):
        self.N, self.M, self.lam = N, M, lam
        self.delay, self.delay_kind = delay, delay_kind
        self.n_blocks = n_blocks
        self.rng = random.Random(seed)
        self.pending_cap = pending_cap
        self.sample_every = sample_every
        # global block store (id 0 = genesis)
        self.lane = [0]
        self.parent = [-1]
        self.link = [-1]
        self.slot = [0]
        self.lh = [0]        # height in its lane chain (genesis = 0)
        self.tmined = [0.0]
        self.miner = [-1]
        self.snap_tips = [None]
        self.snap_c = [0]    # slot of miner's current parent after accepting the block
        self.snap_size = [0]
        self.kdone = [0]     # lanes already holding a child of c when mined (c-children only, else -1)
        # per-node state
        self.slotmap = [{0: 0} for _ in range(M)]
        self.sorted_slots = [[0] for _ in range(M)]
        self.tips = [tuple([0] * N) for _ in range(M)]
        self.cslot = [0] * M
        self.size = [0] * M
        self.pending = [[] for _ in range(M)]  # (recv_time, block)
        self.lane_of_miner = [i % N + 1 for i in range(M)]
        # stats
        self.rej = {'window': 0, 'slot': 0, 'link': 0}
        self._rej_dummy = {'window': 0, 'slot': 0, 'link': 0}
        self.deliveries = 0
        self.adoptions = 0
        self.samples = []   # (time, tuple of node tips)
        self.pending_keep = max(20.0 * delay, 2.0 / lam * N) if delay > 0 else 1.0

    # ---------------- per-node helpers ----------------
    def next_after(self, v, s):
        ss = self.sorted_slots[v]
        i = bisect_right(ss, s)
        return ss[i] if i < len(ss) else None

    def n_children(self, v, s):
        sm = self.slotmap[v]
        base = self.N * s
        return sum(1 for a in range(1, self.N + 1) if (base + a) in sm)

    def try_accept(self, v, b, count=False):
        """Algorithm 1 at node v.  Returns True if accepted."""
        N = self.N
        rej = self.rej if count else self._rej_dummy
        sm = self.slotmap[v]
        p = self.parent[b]
        ps = self.slot[p]
        c = self.cslot[v]
        if sm.get(ps) != p:
            rej['window'] += 1
            return False
        if ps != c:
            n = self.next_after(v, c)
            if ps != n:
                rej['window'] += 1
                return False
        s = self.slot[b]
        if s in sm:
            rej['slot'] += 1
            return False
        a = self.lane[b]
        tips = self.tips[v]
        if self.link[b] != tips[a - 1]:
            rej['link'] += 1
            return False
        # accept
        sm[s] = b
        insort(self.sorted_slots[v], s)
        t = list(tips)
        t[a - 1] = b
        self.tips[v] = tuple(t)
        self.size[v] += 1
        if ps != c:          # parent == n  -> rule 2
            self.cslot[v] = ps
        elif self.n_children(v, c) == N:
            self.cslot[v] = self.next_after(v, c)
        return True

    def adopt(self, v, b):
        """Replace node v's tree by the tree carried in block b (miner's view)."""
        self.adoptions += 1
        new_tips = self.snap_tips[b]
        old_tips = self.tips[v]
        sm = self.slotmap[v]
        ss = self.sorted_slots[v]
        removed, added = [], []
        lh, link = self.lh, self.link
        for a in range(self.N):
            x, y = old_tips[a], new_tips[a]
            while x != y:
                if lh[x] >= lh[y]:
                    removed.append(x)
                    x = link[x]
                else:
                    added.append(y)
                    y = link[y]
        for x in removed:
            s = self.slot[x]
            if sm.get(s) == x:
                del sm[s]
                del ss[bisect_left(ss, s)]
        for y in added:
            s = self.slot[y]
            assert s not in sm
            sm[s] = y
            insort(ss, s)
        self.tips[v] = new_tips
        self.cslot[v] = self.snap_c[b]
        self.size[v] = self.snap_size[b]
        # displaced blocks become candidates for re-application
        if removed:
            pend = self.pending[v]
            now = self.now
            for x in sorted(removed):
                pend.append((now, x))

    def retry_pending(self, v):
        pend = self.pending[v]
        if not pend:
            return
        now = self.now
        keep = self.pending_keep
        sm = self.slotmap[v]
        pend[:] = [(t, x) for (t, x) in pend
                   if now - t <= keep and sm.get(self.slot[x]) != x][-self.pending_cap:]
        progress = True
        while progress and pend:
            progress = False
            rest = []
            for (t, x) in pend:
                if sm.get(self.slot[x]) == x:
                    continue
                if self.try_accept(v, x):
                    progress = True
                else:
                    rest.append((t, x))
            pend[:] = rest

    def deliver(self, v, b):
        self.deliveries += 1
        sm = self.slotmap[v]
        changed = False
        if sm.get(self.slot[b]) != b:
            if self.try_accept(v, b, count=True):
                changed = True
            else:
                self.pending[v].append((self.now, b))
        if self.snap_size[b] > self.size[v]:
            self.adopt(v, b)
            changed = True
        if changed:
            self.retry_pending(v)

    # ---------------- mining ----------------
    def mine(self, u):
        N = self.N
        a = self.lane_of_miner[u]
        sm = self.slotmap[u]
        c = self.cslot[u]
        s = N * c + a
        kd = -1
        if s not in sm:
            p_slot = c
            kd = self.n_children(u, c)
        else:
            n = self.next_after(u, c)
            assert n is not None
            p_slot = n
            s = N * n + a
            assert s not in sm
        b = len(self.lane)
        tip = self.tips[u][a - 1]
        self.lane.append(a)
        self.parent.append(sm[p_slot])
        self.link.append(tip)
        self.slot.append(s)
        self.lh.append(self.lh[tip] + 1)
        self.tmined.append(self.now)
        self.miner.append(u)
        self.kdone.append(kd)
        self.snap_tips.append(None)
        self.snap_c.append(0)
        self.snap_size.append(0)
        ok = self.try_accept(u, b)
        assert ok
        self.snap_tips[b] = self.tips[u]
        self.snap_c[b] = self.cslot[u]
        self.snap_size[b] = self.size[u]
        self.retry_pending(u)
        return b

    def run(self):
        rng = self.rng
        heap = []
        seq = 0
        t = 0.0
        M = self.M
        d = self.delay
        expo = self.delay_kind == 'exp'
        for i in range(self.n_blocks):
            t += rng.expovariate(self.lam)
            while heap and heap[0][0] <= t:
                tt, _, v, b = heapq.heappop(heap)
                self.now = tt
                self.deliver(v, b)
            self.now = t
            if i % self.sample_every == 0:
                # state just before a Poisson mining instant = time-average state (PASTA)
                self.samples.append((t, tuple(self.tips)))
            u = rng.randrange(M)
            b = self.mine(u)
            if d == 0:
                for v in range(M):
                    if v != u:
                        self.deliver(v, b)
            else:
                for v in range(M):
                    if v != u:
                        dd = rng.expovariate(1.0 / d) if expo else d
                        seq += 1
                        heapq.heappush(heap, (t + dd, seq, v, b))
        while heap:
            tt, _, v, b = heapq.heappop(heap)
            self.now = tt
            self.deliver(v, b)
        self.t_end = t
        return self

    # ---------------- analysis ----------------
    def analyze(self, lo_frac=0.05, hi_frac=0.95):
        N = self.N
        # final canonical tree: largest tree; ties -> most common among nodes
        best = max(self.size)
        cands = [v for v in range(self.M) if self.size[v] == best]
        from collections import Counter
        cnt = Counter(self.tips[v] for v in cands)
        ftips = cnt.most_common(1)[0][0]
        fv = next(v for v in cands if self.tips[v] == ftips)
        F = self.slotmap[fv]
        Fids = set(F.values())
        fc = self.cslot[fv]
        nb = len(self.lane) - 1
        lo, hi = int(1 + lo_frac * nb), int(1 + hi_frac * nb)
        t_lo, t_hi = self.tmined[lo], self.tmined[hi]
        mined = hi - lo
        late = race = lanec = casc = 0
        succ = {}
        for x in Fids:
            if x:
                succ[self.link[x]] = x
        for b in range(lo, hi):
            if b in Fids:
                continue
            p, l = self.parent[b], self.link[b]
            if p not in Fids or l not in Fids:
                casc += 1          # built on a block that was itself wasted
                continue
            s = self.slot[b]
            if s in F:
                race += 1          # same-slot race lost (first-seen / fork choice)
                continue
            y = succ.get(l)
            if y is None or self.slot[self.parent[y]] > self.slot[p]:
                late += 1          # its parent was retired before it was accepted (window rule)
            else:
                lanec += 1         # miner's window was ahead of the canonical one ("early")
        wasted = late + race + lanec + casc
        acc = mined - wasted
        # fill rate over retired parents mined inside the window
        par = 0
        ch = 0
        for b in range(lo, hi):
            if b in Fids and self.slot[b] < fc:
                par += 1
                base = N * self.slot[b]
                ch += sum(1 for a in range(1, N + 1) if (base + a) in F)
        fill = ch / (N * par) if par else float('nan')
        # view disagreement at sampled instants
        nd = nfork = 0
        frac_fork = 0.0
        frac_off = 0.0
        ns = 0
        tipcache = {}
        for (t, tl) in self.samples:
            if t < t_lo or t > t_hi:
                continue
            ns += 1
            if any(x != tl[0] for x in tl):
                nd += 1
            from collections import Counter as C2
            mc = C2(tl).most_common(1)[0][1]
            frac_off += 1 - mc / len(tl)
            k = 0
            for x in tl:
                r = tipcache.get(x)
                if r is None:
                    r = any(y not in Fids for y in x)
                    tipcache[x] = r
                k += r
            if k:
                nfork += 1
            frac_fork += k / len(tl)
        # analytic ingredient: mean k over c-children (as mined)
        kc = [self.kdone[b] for b in range(lo, hi)]
        Kbar = sum(max(k, 0) for k in kc) / len(kc)
        return dict(
            mined=mined, wasted=wasted, late=late, race=race, lanec=lanec, casc=casc,
            wasted_frac=wasted / mined, late_frac=(late + lanec) / mined,
            race_frac=race / mined, casc_frac=casc / mined, early_frac=lanec / mined,
            fill=fill, eff_rate=acc / ((t_hi - t_lo) * self.lam),
            p_disagree=nd / ns if ns else float('nan'),
            frac_nodes_offplurality=frac_off / ns if ns else float('nan'),
            p_any_node_on_losing_fork=nfork / ns if ns else float('nan'),
            frac_nodes_on_losing_fork=frac_fork / ns if ns else float('nan'),
            adoptions_per_block=self.adoptions / nb,
            rej_window=self.rej['window'] / self.deliveries,
            rej_slot=self.rej['slot'] / self.deliveries,
            rej_link=self.rej['link'] / self.deliveries,
            Kbar=Kbar,
        )
