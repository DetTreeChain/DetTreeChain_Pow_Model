import time, sys
from dtsim import Sim
for N in [1,2,3,4,8]:
    for d in [0.0, 0.1]:
        t0=time.time()
        s = Sim(N, 4*N, N/1.0, d, 'fixed', 10000, 1).run()
        r = s.analyze()
        print(N, d, f"{time.time()-t0:.1f}s", {k: round(v,4) if isinstance(v,float) else v for k,v in r.items()})
