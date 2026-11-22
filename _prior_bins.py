import csv, numpy as np
rows=[r for r in csv.DictReader(open("data/alcdef_identifiability_catalog.csv",encoding="utf-8"))]
def f(r,k):
    try: return float(r[k])
    except: return float("nan")
main=[r for r in rows if np.isfinite(f(r,"r4")) and np.isfinite(f(r,"rstar")) and f(r,"A1")>=0.015 and f(r,"A1sig")>=1.5]
r4=np.array([f(r,"r4") for r in main]); rs=np.array([f(r,"rstar") for r in main])
d2=np.array([f(r,"d2") for r in main]); d4=np.array([f(r,"d4") for r in main])
danger=(r4>rs)&((d2<0)&(d4>0))
print("n =", len(main))
for t0 in (0.15,0.20,0.25,0.80):
    pred=100*np.mean(r4>t0)
    m=(rs>=t0-0.025)&(rs<t0+0.025)
    obs=100*np.mean(danger[m]) if m.sum() else float("nan")
    print("r0*=%.2f  predicted %5.1f%%  observed %5.1f%%  n=%d" % (t0,pred,obs,m.sum()))
