import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import analytic

rows = list(csv.DictReader(open('results.csv')))
for r in rows:
    for k, v in r.items():
        try:
            r[k] = float(v)
        except (ValueError, TypeError):
            pass

COL = {1: '#2a78d6', 2: '#eb6834', 3: '#1baf7a', 4: '#eda100', 8: '#4a3aa7'}
NS = [1, 2, 3, 4, 8]
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#888', 'axes.linewidth': 0.6,
                     'xtick.color': '#52514e', 'ytick.color': '#52514e'})
fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), facecolor='#fcfcfb')
xs_fine = np.linspace(0, 0.2, 81)


def get(design, N, kind, m='wasted_frac'):
    rs = sorted([r for r in rows if r['design'] == design and r['N'] == N and r['kind'] == kind],
                key=lambda r: r['x'])
    return (np.array([r['x'] for r in rs]), np.array([r[m] for r in rs]),
            np.array([r[m + '_se'] for r in rs]))


for j, kind in enumerate(['fixed', 'exp']):
    ax = axes[0, j]
    for N in NS:
        x, y, e = get('DTC', N, kind)
        ax.errorbar(x, y, yerr=2 * e, color=COL[N], lw=2, marker='o', ms=4, capsize=2, label=f'N={N} sim')
        ya = [analytic.w_dtc(N, t)[2] for t in xs_fine]
        ax.plot(xs_fine, ya, color=COL[N], lw=1, ls=':')
        ax.annotate(f'N={N}', (x[-1], y[-1]), xytext=(4, 0), textcoords='offset points',
                    color='#0b0b0b', va='center', fontsize=8)
    ax.set_title(f'DetTreeChain wasted fraction ({kind} delay)\nsolid = simulation (±2 SE), dotted = analytic approx.',
                 fontsize=9.5, loc='left')
    ax.set_xlabel('Δ / I   (I = per-lane block interval)')
    ax.set_ylabel('fraction of mined blocks not in final tree')
    ax.set_xlim(0, 0.225)
    ax.grid(color='#e5e5e2', lw=0.5)
    ax.legend(fontsize=7.5, frameon=False, loc='upper left')

    ax = axes[1, j]
    for N in NS:
        x, y, e = get('DTC', N, kind)
        ax.plot(x[1:], y[1:], color=COL[N], lw=2, marker='o', ms=4, label=f'N={N} DetTreeChain')
        if N > 1:
            x1, y1, _ = get('B1', N, kind)
            ax.plot(x1[1:], y1[1:], color=COL[N], lw=1.5, ls='--', marker='s', ms=3)
            x2, y2, _ = get('B2', N, kind)
            ax.plot(x2[1:], y2[1:], color=COL[N], lw=1.2, ls=':', marker='^', ms=3)
    ax.plot([], [], color='#52514e', lw=2, label='solid: DetTreeChain (N lanes, rate N/I)')
    ax.plot([], [], color='#52514e', lw=1.5, ls='--', label='dashed: B1 chain, interval I/N, delay Δ')
    ax.plot([], [], color='#52514e', lw=1.2, ls=':', label='dotted: B2 chain, interval I, delay NΔ (≈B1)')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_title(f'vs. single-chain baselines at the same aggregate rate ({kind} delay)\ncolor = N (N=1: all three coincide)',
                 fontsize=9.5, loc='left')
    ax.set_xlabel('Δ / I')
    ax.set_ylabel('wasted / stale fraction (log)')
    ax.grid(color='#e5e5e2', lw=0.5, which='both')
    ax.legend(fontsize=7, frameon=False, loc='lower right', ncol=1)
for ax in axes.flat:
    ax.set_facecolor('#fcfcfb')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig('wasted_vs_delay.png', dpi=150)
print('saved')
