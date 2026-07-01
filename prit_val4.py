"""
4-point model extrapolation validation.
Uses the 4-point power law fit (trained on 10^6-10^9) to predict pi_2(10^10) and pi_2(10^11).
"""
import math
import numpy as np
from scipy.optimize import minimize

C2 = 0.6601618914

pi2 = {
    10**6:  8169,
    10**7:  58980,
    10**8:  440312,
    10**9:  3424506,
    10**10: 27412679,
    10**11: 224376048,
}

def hl(x):
    return 2 * C2 * x / math.log(x)**2

def residual(x):
    return pi2[x] - hl(x)

# Data
x4 = np.array([1e6, 1e7, 1e8, 1e9])
r4 = np.array([residual(int(x)) for x in x4])

# Simple power law fit: R(x) = C * x^alpha
log_x = np.log10(x4)
log_R = np.log10(r4)
coeffs = np.polyfit(log_x, log_R, 1)
alpha_4 = coeffs[0]
C_4 = 10**coeffs[1]

print("=== 4-POINT POWER LAW FIT ===")
print(f"R(x) = {C_4:.4e} * x^{alpha_4:.4f}")
print()

# Evaluate on training data
ss_res = 0
for i, x in enumerate(x4):
    pred = C_4 * x**alpha_4
    err = abs(pred - r4[i]) / r4[i] * 100
    ss_res += (pred - r4[i])**2
    print(f"  x={int(x):>10d}  R={r4[i]:>12.1f}  pred={pred:>12.1f}  err={err:>6.1f}%")

r2 = 1 - ss_res / sum((r4 - np.mean(r4))**2)
print(f"\nR^2 = {r2:.10f}")
print(f"RMSE = {math.sqrt(ss_res):.1f}")

# Extrapolate to 10^10
print("\n=== EXTRAPOLATION TO 10^10 ===")
x_pred = 10**10
R_pred = C_4 * x_pred**alpha_4
R_act = residual(x_pred)
hl_val = hl(x_pred)
pi2_pred = hl_val + R_pred
pi2_act = pi2[x_pred]
err_R = abs(R_pred - R_act) / R_act * 100
err_pi2 = abs(pi2_pred - pi2_act) / pi2_act * 100
print(f"  R_pred = {R_pred:.1f}")
print(f"  R_act  = {R_act:.1f}")
print(f"  R_err  = {err_R:.1f}%")
print(f"  pi2_pred = {pi2_pred:,.0f}")
print(f"  pi2_act  = {pi2_act:,d}")
print(f"  pi2_err  = {err_pi2:.2f}%")

# Extrapolate to 10^11
print("\n=== EXTRAPOLATION TO 10^11 ===")
x_pred = 10**11
R_pred = C_4 * x_pred**alpha_4
R_act = residual(x_pred)
hl_val = hl(x_pred)
pi2_pred = hl_val + R_pred
pi2_act = pi2[x_pred]
err_R = abs(R_pred - R_act) / R_act * 100
err_pi2 = abs(pi2_pred - pi2_act) / pi2_act * 100
print(f"  R_pred = {R_pred:.1f}")
print(f"  R_act  = {R_act:.1f}")
print(f"  R_err  = {err_R:.1f}%")
print(f"  pi2_pred = {pi2_pred:,.0f}")
print(f"  pi2_act  = {pi2_act:,d}")
print(f"  pi2_err  = {err_pi2:.2f}%")

# Note
print("\n=== CAVEAT ===")
print("Extrapolation beyond 10^11 assumes the power law continues to hold.")
print("This has not been verified computationally.")
