"""
Power law analysis of the Hardy-Littlewood residual R(x) = pi_2(x) - 2*C2*x/(log x)^2.
Computes growth ratios, local exponents, and fits R(x) = C * x^alpha.
"""
import math
import numpy as np
from scipy.optimize import curve_fit

C2 = 0.6601618914

# Verified twin prime counts
pi2 = {
    10**6:  8169,
    10**7:  58980,
    10**8:  440312,
    10**9:  3424506,
    10**10: 27412679,
    10**11: 224376048,
}

def hl(x):
    """Hardy-Littlewood main term 2*C2*x/(log x)^2."""
    return 2 * C2 * x / math.log(x)**2

def residual(x):
    """R(x) = pi_2(x) - HL(x)."""
    return pi2[x] - hl(x)

# Compute residuals
print("=== RESIDUALS ===")
print(f"{'x':>12s}  {'pi_2(x)':>12s}  {'HL(x)':>14s}  {'R(x)':>14s}")
xs = sorted(pi2.keys())
for x in xs:
    r = residual(x)
    print(f"  {x:>10d}  {pi2[x]:>12d}  {hl(x):>14.1f}  {r:>14.1f}")

# Growth ratios
print("\n=== GROWTH RATIOS ===")
print(f"{'x_i':>12s}  {'x_{i+1}':>12s}  {'R(x_i)':>12s}  {'R(x_{i+1})':>12s}  {'ratio':>8s}  {'alpha_loc':>10s}")
R = [residual(x) for x in xs]
for i in range(len(xs) - 1):
    ratio = R[i+1] / R[i]
    alpha_loc = math.log10(ratio)
    print(f"  {xs[i]:>10d}  {xs[i+1]:>10d}  {R[i]:>12.1f}  {R[i+1]:>12.1f}  {ratio:>8.2f}  {alpha_loc:>10.3f}")

# Simple power law fit: log(R) = log(C) + alpha * log(x)
print("\n=== LOG-LOG LINEAR FIT ===")
log_x = np.log10(xs)
log_R = np.log10(R)

# Fit all 6 points
coeffs = np.polyfit(log_x, log_R, 1)
alpha_all = coeffs[0]
C_all = 10**coeffs[1]
print(f"6-point fit: R(x) = {C_all:.4e} * x^{alpha_all:.4f}")

# Fit first 4 points
coeffs4 = np.polyfit(log_x[:4], log_R[:4], 1)
alpha_4 = coeffs4[0]
C_4 = 10**coeffs4[1]
print(f"4-point fit: R(x) = {C_4:.4e} * x^{alpha_4:.4f}")

# Fit first 5 points
coeffs5 = np.polyfit(log_x[:5], log_R[:5], 1)
alpha_5 = coeffs5[0]
C_5 = 10**coeffs5[1]
print(f"5-point fit: R(x) = {C_5:.4e} * x^{alpha_5:.4f}")

# Evaluate fits
def power_law(x, C, alpha):
    return C * x**alpha

print("\n=== FIT EVALUATION ===")
for name, C_f, a_f, n_pts in [("6-point", C_all, alpha_all, 6), ("4-point", C_4, alpha_4, 4), ("5-point", C_5, alpha_5, 5)]:
    print(f"\n{name} fit (trained on {n_pts} points):")
    ss_res = 0
    ss_tot = 0
    for i, x in enumerate(xs):
        pred = power_law(x, C_f, a_f)
        err = abs(pred - R[i]) / R[i] * 100
        ss_res += (pred - R[i])**2
        ss_tot += (R[i] - np.mean(R))**2
        marker = " *" if i < n_pts else ""
        print(f"  x={x:>10d}  R={R[i]:>12.1f}  pred={pred:>12.1f}  err={err:>6.1f}%{marker}")
    r2 = 1 - ss_res / ss_tot
    rmse = math.sqrt(ss_res)
    print(f"  R^2 = {r2:.8f}, RMSE = {rmse:.1f}")

# Extrapolation
print("\n=== EXTRAPOLATION ===")
for name, C_f, a_f in [("4-point", C_4, alpha_4), ("5-point", C_5, alpha_5)]:
    for xp in [10**10, 10**11]:
        R_pred = power_law(xp, C_f, a_f)
        R_act = residual(xp)
        hl_val = hl(xp)
        pi2_pred = hl_val + R_pred
        pi2_act = pi2[xp]
        err_R = abs(R_pred - R_act) / R_act * 100
        err_pi2 = abs(pi2_pred - pi2_act) / pi2_act * 100
        print(f"  {name} at x=10^{int(math.log10(xp))}:")
        print(f"    R_pred={R_pred:>14.1f}  R_act={R_act:>14.1f}  R_err={err_R:.1f}%")
        print(f"    pi2_pred={pi2_pred:>14.0f}  pi2_act={pi2_act:>14d}  pi2_err={err_pi2:.2f}%")

# R^2 for the 5-point fit
print("\n=== R^2 FOR 5-POINT FIT ===")
ss_res = sum((power_law(x, C_5, alpha_5) - residual(x))**2 for x in xs[:5])
ss_tot = sum((residual(x) - np.mean([residual(x) for x in xs[:5]]))**2 for x in xs[:5])
print(f"R^2 = {1 - ss_res/ss_tot:.10f}")
