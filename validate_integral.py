"""
VALIDATION: Follow the other AI's 3-step method exactly.
Step 1: R_int(x) = pi_2(x) - 2C2 * (Li(x) - x/log(x))
Step 2: Refit power law
Step 3: Compare
"""
import math
from scipy.integrate import quad
from scipy.optimize import curve_fit
import numpy as np

C2 = 0.6601618914

nicely_data = {
    10**6: 8169, 10**7: 58980, 10**8: 440312, 10**9: 3424506,
    10**10: 27412679, 10**11: 224376048,
    2e11: 424084653, 3e11: 615885700, 4e11: 802817718, 5e11: 986222314,
    6e11: 1166916933, 7e11: 1345394380, 8e11: 1521998439, 9e11: 1696987738,
    10**12: 1870585220,
    2e12: 3552770943, 3e12: 5173760785, 4e12: 6756832076, 5e12: 8312493003,
    6e12: 9846842484, 7e12: 11363874338, 8e12: 12866256870, 9e12: 14356002120,
    10**13: 15834664872,
    2e13: 30198862775, 3e13: 44078684643, 4e13: 57657248284, 5e13: 71018282471,
    6e13: 84209699420, 7e13: 97262712867, 8e13: 110198743491, 9e13: 123033833767,
    10**14: 135780321665,
}

def li(x):
    """Logarithmic integral Li(x) = integral from 2 to x of dt/log(t)"""
    result, _ = quad(lambda t: 1/math.log(t), 2, x)
    return result

def hl_simplified(x):
    """2C2 * x / (log x)^2"""
    return 2 * C2 * x / math.log(x)**2

def hl_integral(x):
    """2C2 * (Li(x) - x/log(x))"""
    return 2 * C2 * (li(x) - x/math.log(x))

def power_law(x, C, alpha):
    return C * x**alpha

print("=" * 72)
print("STEP 1: Recompute R(x) with integral")
print("=" * 72)
print(f"R_int(x) = pi_2(x) - 2C2 * (Li(x) - x/log(x))")
print()
print(f"{'x':>8} {'pi_2(x)':>15} {'HL_simp':>15} {'HL_int':>15} {'R_simp':>12} {'R_int':>12}")
print("-" * 82)

xs = sorted(nicely_data.keys())
R_simp_list = []
R_int_list = []

for x in xs:
    pi2 = nicely_data[x]
    hs = hl_simplified(x)
    hi = hl_integral(x)
    rs = pi2 - hs
    ri = pi2 - hi
    R_simp_list.append(rs)
    R_int_list.append(ri)
    print(f"{x:8.0e} {pi2:15d} {hs:15.1f} {hi:15.1f} {rs:12.1f} {ri:12.1f}")

xs = np.array(xs)
R_simp = np.array(R_simp_list)
R_int = np.array(R_int_list)

print(f"\nR_simp range: [{R_simp.min():.1f}, {R_simp.max():.1f}]")
print(f"R_int range:  [{R_int.min():.1f}, {R_int.max():.1f}]")
print(f"R_int/R_simp ratio: {np.mean(np.abs(R_int)/np.abs(R_simp)):.4f}")

print("\n" + "=" * 72)
print("STEP 2: Refit power law to R_int")
print("=" * 72)

# Fit simplified (for comparison)
popt_simp, _ = curve_fit(power_law, xs, R_simp, p0=[0.01, 0.86], maxfev=10000)
pred_simp = power_law(xs, *popt_simp)
r2_simp = 1 - np.sum((R_simp - pred_simp)**2)/np.sum((R_simp - np.mean(R_simp))**2)

# Fit integral
popt_int, _ = curve_fit(power_law, xs, R_int, p0=[0.01, 0.86], maxfev=10000)
pred_int = power_law(xs, *popt_int)
r2_int = 1 - np.sum((R_int - pred_int)**2)/np.sum((R_int - np.mean(R_int))**2)

print(f"\nSimplified: R(x) = {popt_simp[0]:.4e} * x^{popt_simp[1]:.4f}")
print(f"  R^2 = {r2_simp:.6f}")
print(f"  alpha = {popt_simp[1]:.4f}")

print(f"\nIntegral:   R(x) = {popt_int[0]:.4e} * x^{popt_int[1]:.4f}")
print(f"  R^2 = {r2_int:.6f}")
print(f"  alpha = {popt_int[1]:.4f}")

print("\n" + "=" * 72)
print("STEP 3: Compare")
print("=" * 72)

print(f"\nAlpha difference: {abs(popt_simp[1] - popt_int[1]):.4f}")
print(f"R^2 difference:   {abs(r2_simp - r2_int):.4f}")

if r2_int > 0.95:
    print("\n>>> OUTCOME: Power law SURVIVES (alpha ~ 0.86)")
    print(">>> Your discovery is ROBUST")
elif r2_int > 0.5:
    print("\n>>> OUTCOME: Power law WEAKENS (alpha ~ 0.5-0.6)")
    print(">>> Your result was PARTIALLY an artifact")
else:
    print("\n>>> OUTCOME: Power law DISAPPEARS")
    print(">>> Your result was ENTIRELY an artifact")
    print(f">>> R^2 dropped from {r2_simp:.4f} to {r2_int:.4f}")

# Also check: what does the integral residual look like?
print("\n" + "=" * 72)
print("SIGN ANALYSIS")
print("=" * 72)
signs_simp = np.sign(R_simp)
signs_int = np.sign(R_int)
print(f"Simplified: all positive = {np.all(signs_simp > 0)}")
print(f"Integral:   all positive = {np.all(signs_int > 0)}")
print(f"Integral signs: {np.sum(signs_int > 0)} positive, {np.sum(signs_int < 0)} negative")
