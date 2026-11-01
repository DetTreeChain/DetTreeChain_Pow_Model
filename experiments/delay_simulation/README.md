# DetTreeChain: delay simulation (Section 8.4 / RQ3 of the paper)

Discrete-event simulation of the DetTreeChain window rule (Rules 1 and 2) under propagation delay,
compared with a single chain at the same total block rate.

## Model
- N lanes, 4 miners per lane, equal hash power, total block rate N/I (I = interval of one lane, set to 1).
- Full mesh; every (block, node) delay is either fixed (Delta) or exponential with mean Delta.
- Nodes apply the acceptance rules in arrival order (first seen wins) and adopt a peer's tree when it holds
  strictly more blocks (a simplified longest-tree rule).
- A block counts as wasted if it is not in the final tree.
- Baselines (same code with N = 1): B1 = one chain with interval I/N; B2 = one chain with interval I and
  N-times larger blocks, modelled as delay N*Delta.

## Files
| File | Purpose |
| --- | --- |
| `dtsim.py` | The simulator (class `Sim`). |
| `analytic.py` | Analytic approximation and the zero-delay window model (fill rate, lanes under c'). |
| `run_all.py` | Full sweep: N in {1,2,3,4,8}, Delta/I in {0,0.01,0.02,0.05,0.1,0.2}, fixed and exponential delay, 5 seeds x 50,000 blocks. Writes `results_raw.csv` (per seed) and `results.csv` (mean and standard error). |
| `quick.py` | Short sanity run. |
| `plot.py` | Draws `wasted_vs_delay.png` from `results.csv`. |
| `results.csv`, `results_raw.csv` | The results reported in the paper. |

## Running
```
python3 run_all.py          # full sweep, uses all CPU cores (about 1-2 minutes on 20+ cores)
pip install -r requirements.txt && python3 plot.py
```
The paper's Figure (wasted blocks vs. Delta/I) uses the `fixed` rows of `results.csv`
(column `wasted_frac`) for designs `DTC` and `B1`.
