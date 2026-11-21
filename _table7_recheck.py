import csv, numpy as np
rows=[r for r in csv.DictReader(open("data/alcdef_identifiability_catalog.csv",encoding="utf-8"))]
def f(r,k):
    try: return float(r[k])
    except: return float("nan")
main=[r for r in rows if np.isfinite(f(r,"r4")) and np.isfinite(f(r,"rstar")) and f(r,"A1")>=0.015 and f(r,"A1sig")>=1.5]
r4=np.array([f(r,"r4") for r in main]); rs=np.array([f(r,"rstar") for r in main])
danger = r4 > rs
print(" t0   pred   obs(fig method)  n")
for t0 in (0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.70,0.80,0.90):
    m=(rs>=t0-0.025)&(rs<t0+0.025)
    print("%.2f  %5.1f   %5.1f            %d" % (t0, 100*np.mean(r4>t0), 100*np.mean(danger[m]) if m.sum() else float("nan"), m.sum()))
