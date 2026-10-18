"""Analytic approximation for DetTreeChain wasted work under small delay.

Zero-delay window process (equal lanes, each lane rate 1/I): a cycle of the current
parent c starts with k0 lanes already holding a child of c (k0=1 if the previous parent
was retired by Rule 2, k0=0 if it was retired full).  While k lanes are done, the next
block is a child of c w.p. (N-k)/N (k -> k+1) and the first child of n w.p. k/N (cycle ends).

A child of c mined while k lanes are done is *late* if a first child of n was mined in
the preceding delay window and had not reached its miner: P ~ 1 - exp(-k*Delta/I).
With fixed delay every other node then rejects it; it survives only if its own miner
(1 of m miners of the lane) finds one of the next two lane blocks first, so it is lost
w.p. ((m-1)/m)^2.  Same-slot races (two miners of one lane within Delta) cost
(m-1)/m * (1 - exp(-Delta/I)) per block, exactly as in a chain of rate 1/I.

  w_DTC(Delta) ~ (m-1)/m * (1-e^{-x}) + ((m-1)/m)^2 * sum_k pi_k (1 - e^{-k x}),   x = Delta/I
               ~ [ (m-1)/m + ((m-1)/m)^2 * K_N ] * x                (first order)
  K_N = E[k * 1{block is a child of c}]  (per mined block)

Baselines (chain, M miners): B1 rate N/I, delay Delta; B2 rate 1/I, delay N*Delta:
  w_B1 ~ w_B2 ~ (M-1)/M * (1 - e^{-N x})  ~ N x.
"""
import math


def window_chain(N):
    """Return (fill, pi) with pi[k] = P(a mined block is a child of c mined with k lanes done)."""
    if N == 1:
        return 1.0, {}
    # stationary distribution of k0 in {0,1}
    def cycle(k0):
        # returns P(end by rule 2), expected c-children, expected per-k counts
        probs = {k0: 1.0}
        endr2 = 0.0
        cch = 0.0
        perk = [0.0] * N
        full = 0.0
        k = k0
        p = 1.0
        # path is deterministic in k: from state k either c-child (k+1) or end
        while True:
            if k == N:
                full = p
                break
            pc = (N - k) / N
            pn = k / N
            perk[k] += p * pc
            cch += p * pc
            endr2 += p * pn
            p *= pc
            k += 1
        return endr2, full, cch, perk
    e0 = cycle(0); e1 = cycle(1)
    # transition k0: from a cycle, next k0 = 1 w.p. endr2, 0 w.p. full
    # P(k0=1) = q solves q = q*e1[0] + (1-q)*e0[0]
    q = e0[0] / (1 - e1[0] + e0[0])
    blocks = (1 - q) * (e0[2] + e0[0]) + q * (e1[2] + e1[0])   # c-children + ending n-child
    children = (1 - q) * e0[2] + q * (1 + e1[2])
    fill = children / N
    pi = {k: ((1 - q) * e0[3][k] + q * e1[3][k]) / blocks for k in range(N)}
    return fill, pi


def K(N):
    _, pi = window_chain(N)
    return sum(k * p for k, p in pi.items())


def w_dtc(N, x, m=4):
    _, pi = window_chain(N)
    a = (m - 1) / m
    race = a * (1 - math.exp(-x))
    late = a * a * sum(p * (1 - math.exp(-k * x)) for k, p in pi.items())
    return race, late, race + late


def w_chain(rate_x, M):
    """stale fraction of a chain with M equal miners, rate*delay = rate_x"""
    return (M - 1) / M * (1 - math.exp(-rate_x))


if __name__ == '__main__':
    for N in [1, 2, 3, 4, 8, 16, 32, 64, 128]:
        f, pi = window_chain(N)
        print(N, 'fill=%.4f' % f, 'K_N=%.4f' % K(N), 'K/sqrtN=%.3f' % (K(N) / math.sqrt(N)))
