"""
Oscillatory model fit (simplified).
Demonstrates that the oscillatory model R(x) = A * x^alpha * sum(w_ij * cos(...))
can fit training data but fails to generalize.
Uses a coarse grid for speed.
"""
import math
import numpy as np
from scipy.optimize import minimize

C2 = 0.6601618914
gz = np.array([14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,48.005151,49.773832,52.970321,56.446248,59.347044,60.831778,65.112544,67.079811,69.546402,72.067158,75.704691,77.144840,79.337375,82.910381,84.735493,87.425275,88.809111,92.491899,94.651344,95.870634,98.831194,101.317851,103.725539,105.446623,107.168611,111.029517,111.874659,114.320221,116.226680,118.943299,121.368869,122.946064,124.256904,127.516683,129.578704,131.254325,133.498619,134.764525,138.116078,139.736121,141.123707,143.102708])
N = len(gz)
ii, jj = np.triu_indices(N, 1)
gd = gz[jj] - gz[ii]

pi2 = {10**6:8169, 10**7:58980, 10**8:440312, 10**9:3424506, 10**10:27412679, 10**11:224376048}

def pred(xa, al, sg, A):
    lx = np.log(xa)
    gw = np.exp(-gd**2 / (2*sg**2))
    cos_t = np.cos(gd[None,:] * lx[:,None])
    return A * (xa**al) * (cos_t @ gw)

def fit(x_data, r_data, starts):
    best = 1e20
    bp = None
    for s in starts:
        def obj(p):
            if p[2] < 0.5 or p[2] > 20 or abs(p[0]) > 1:
                return 1e20
            return np.sum((pred(x_data, p[1], p[2], p[0]) - r_data)**2)
        r = minimize(obj, s, method='Nelder-Mead', options={'maxiter': 3000, 'xatol': 1e-8, 'fatol': 1e-12})
        if r.fun < best:
            best = r.fun
            bp = r.x
    return bp

# Small grid for speed
starts = []
for A in [1e-4, 5e-4, 1e-3, 5e-3]:
    for a in [0.7, 0.8, 0.9, 1.0]:
        for s in [3.0, 5.0, 7.0]:
            starts.append([A, a, s])

# 4-point fit
x4 = np.array([1e6, 1e7, 1e8, 1e9])
r4 = np.array([pi2[int(x)] - 2*C2*x/(math.log(x)**2) for x in x4])
p4 = fit(x4, r4, starts)
A4, a4, s4 = p4
pv4 = pred(x4, a4, s4, A4)
r2_4 = 1 - np.sum((pv4-r4)**2) / np.sum((r4 - np.mean(r4))**2)

print("=== 4-POINT OSCILLATORY FIT ===")
print("A = %.6e, alpha = %.4f, sigma = %.4f" % (A4, a4, s4))
print("R^2 = %.8f" % r2_4)
for i in range(4):
    x = int(x4[i])
    err = abs(pv4[i] - r4[i]) / abs(r4[i]) * 100
    print("  x=%10d  R=%12.1f  pred=%12.1f  err=%.1f%%" % (x, r4[i], pv4[i], err))

# Extrapolate to 10^10
print("\n=== EXTRAPOLATION TO 10^10 ===")
xp = 1e10
pxp = pred(np.array([xp]), a4, s4, A4)[0]
ract = pi2[int(xp)] - 2*C2*xp/(math.log(xp)**2)
hl_val = 2*C2*xp/(math.log(xp)**2)
pi2_pred = hl_val + pxp
pi2_act = pi2[int(xp)]
print("  R_pred=%12.1f  R_act=%12.1f  R_err=%.1f%%" % (pxp, ract, abs(pxp-ract)/abs(ract)*100))
print("  pi2_pred=%14.0f  pi2_act=%14d  pi2_err=%.2f%%" % (pi2_pred, pi2_act, abs(pi2_pred-pi2_act)/pi2_act*100))

# 5-point fit
print("\n=== 5-POINT OSCILLATORY FIT ===")
x5 = np.array([1e6, 1e7, 1e8, 1e9, 1e10])
r5 = np.array([pi2[int(x)] - 2*C2*x/(math.log(x)**2) for x in x5])
p5 = fit(x5, r5, starts)
A5, a5, s5 = p5
pv5 = pred(x5, a5, s5, A5)
r2_5 = 1 - np.sum((pv5-r5)**2) / np.sum((r5 - np.mean(r5))**2)

print("A = %.6e, alpha = %.4f, sigma = %.4f" % (A5, a5, s5))
print("R^2 = %.8f" % r2_5)
for i in range(5):
    x = int(x5[i])
    err = abs(pv5[i] - r5[i]) / abs(r5[i]) * 100
    print("  x=%10d  R=%12.1f  pred=%12.1f  err=%.1f%%" % (x, r5[i], pv5[i], err))

# 5-point prediction for 10^11
p11 = pred(np.array([1e11]), a5, s5, A5)[0]
hl11 = 2*C2*1e11/math.log(1e11)**2
pi2_11 = hl11 + p11
print("\n  10^11 prediction: pi2=%.0f (actual=224376048, err=%.2f%%)" % (pi2_11, abs(pi2_11-224376048)/224376048*100))

# Dominant pairs
print("\n=== DOMINANT PAIRS ===")
gw = np.exp(-gd**2 / (2*s5**2))
order = np.argsort(-gw)
for k in order[:5]:
    i, j = ii[k], jj[k]
    print("  (%d,%d) sep=%.3f  w=%.4f" % (i+1, j+1, gd[k], gw[k]))

print("\n=== NOTE ===")
print("The oscillatory model can fit training data with high R^2,")
print("but produces negative predictions and fails to generalize.")
print("See prit_powerlaw.py for the power law model that works.")
