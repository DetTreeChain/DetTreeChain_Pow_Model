"""Run all configurations; write results_raw.csv (per seed) and results.csv (aggregated)."""
import csv, math, sys, time
from multiprocessing import Pool
from dtsim import Sim
import analytic

NS = [1, 2, 3, 4, 8]
XS = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2]          # Delta / I  (I = per-lane interval = 1)
KINDS = ['fixed', 'exp']
SEEDS = list(range(5))
BLOCKS = 50000
MPL = 4                                          # miners per lane

def jobs():
    for kind in KINDS:
        for x in XS:
            for N in NS:
                M = MPL * N
                for s in SEEDS:
                    yield ('DTC', N, x, kind, s, dict(N=N, M=M, lam=float(N), delay=x))
                    if N > 1:
                        # B1: one chain, block interval I/N, same delay, same miners
                        yield ('B1', N, x, kind, s, dict(N=1, M=M, lam=float(N), delay=x))
                        # B2: one chain, interval I, N-times larger blocks: delay N*Delta
                        yield ('B2', N, x, kind, s, dict(N=1, M=M, lam=1.0, delay=N * x))

def work(job):
    design, N, x, kind, seed, p = job
    import zlib
    sd = zlib.crc32(repr((design, N, x, kind, seed)).encode())
    sim = Sim(p['N'], p['M'], p['lam'], p['delay'], kind, BLOCKS, sd).run()
    r = sim.analyze()
    r.update(design=design, N=N, x=x, kind=kind, seed=seed)
    return r

if __name__ == '__main__':
    t0 = time.time()
    J = list(jobs())
    with Pool(22) as pool:
        res = pool.map(work, J, chunksize=1)
    print('runs', len(res), 'time %.0fs' % (time.time() - t0))
    keys = ['design', 'N', 'x', 'kind', 'seed'] + [k for k in res[0] if k not in ('design', 'N', 'x', 'kind', 'seed')]
    with open('results_raw.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(res)
    # aggregate
    from collections import defaultdict
    g = defaultdict(list)
    for r in res:
        g[(r['design'], r['N'], r['x'], r['kind'])].append(r)
    metrics = ['wasted_frac', 'late_frac', 'early_frac', 'race_frac', 'casc_frac', 'fill',
               'eff_rate', 'p_disagree', 'frac_nodes_offplurality', 'p_any_node_on_losing_fork',
               'frac_nodes_on_losing_fork', 'adoptions_per_block', 'rej_window', 'rej_slot', 'rej_link', 'Kbar']
    out = []
    for (d, N, x, kind), rs in sorted(g.items(), key=lambda t: (t[0][3], t[0][0], t[0][1], t[0][2])):
        row = dict(design=d, N=N, x=x, kind=kind, seeds=len(rs), mined=sum(r['mined'] for r in rs))
        for m in metrics:
            v = [r[m] for r in rs]
            mu = sum(v) / len(v)
            se = (sum((a - mu) ** 2 for a in v) / (len(v) - 1)) ** 0.5 / len(v) ** 0.5
            row[m] = mu; row[m + '_se'] = se
        row['acc_per_mined'] = 1 - row['wasted_frac']
        if d == 'DTC':
            ra, la, to = analytic.w_dtc(N, x, MPL)
            row['analytic_race'] = ra; row['analytic_late'] = la; row['analytic_wasted'] = to
        elif d == 'B1':
            row['analytic_wasted'] = analytic.w_chain(N * x, MPL * N)
        else:
            row['analytic_wasted'] = analytic.w_chain(N * x, MPL * N)
        out.append(row)
    keys = list(out[0].keys())
    for r in out:
        for k in r:
            if k not in keys: keys.append(k)
    with open('results.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(out)
