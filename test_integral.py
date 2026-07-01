"""
TEST: Does the power law survive when using the FULL Hardy-Littlewood
integral instead of the simplified 2C2*x/(log x)^2?
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

def hl_simplified(x):
    """2C2 * x / (log x)^2"""
    return 2 * C2 * x / math.log(x)**2

def hl_integral(x):
    """2C2 * integral from 2 to x of dt/(log t)^2"""
    result, _ = quad(lambda t: 1/math.log(t)**2, 2, x)
    return 2 * C2 * result

print("=" * 72)
print("COMPARISON: Simplified vs Integral Hardy-Littlewood")
print("=" * 72)
print(f"{'x':>8} {'R_simplified':>15} {'R_integral':>15} {'ratio':>8}")
print("-" * 50)

for x in sorted(nicely_data.keys()):
    pi2 = nicely_data[x]
    r_simp = pi2 - hl_simplified(x)
    r_int = pi2 - hl_integral(x)
    ratio = r_simp / r_int if r_int != 0 else float('inf')
    print(f"{x:8.0e} {r_simp:15.1f} {r_int:15.1f} {ratio:8.4f}")

# Now fit power law to BOTH residuals
print("\n" + "=" * 72)
print("POWER LAW FITS")
print("=" * 72)

def power_law(x, C, alpha):
    return C * x**alpha

xs = np.array(sorted(nicely_data.keys()))
pi2s = np.array([nicely_data[int(x)] for x in xs])

# Simplified
R_simp = np.array([pi2s[i] - hl_simplified(xs[i]) for i in range(len(xs))])
# Integral
R_int = np.array([pi2s[i] - hl_integral(xs[i]) for i in range(len(xs))])

# Fit simplified
popt_simp, _ = curve_fit(power_law, xs, R_simp, p0=[0.01, 0.86], maxfev=10000)
pred_simp = power_law(xs, *popt_simp)
ss_res_simp = np.sum((R_simp - pred_simp)**2)
ss_tot_simp = np.sum((R_simp - np.mean(R_simp))**2)
r2_simp = 1 - ss_res_simp/ss_tot_simp

# Fit integral
popt_int, _ = curve_fit(power_law, xs, R_int, p0=[0.01, 0.86], maxfev=10000)
pred_int = power_law(xs, *popt_int)
ss_res_int = np.sum((R_int - pred_int)**2)
ss_tot_int = np.sum((R_int - np.mean(R_int))**2)
r2_int = 1 - ss_res_int/ss_tot_int

print(f"\nSimplified: R(x) = {popt_simp[0]:.4e} * x^{popt_simp[1]:.4f}")
print(f"  R^2 = {r2_simp:.6f}")
print(f"  alpha = {popt_simp[1]:.4f}")

print(f"\nIntegral:   R(x) = {popt_int[0]:.4e} * x^{popt_int[1]:.4f}")
print(f"  R^2 = {r2_int:.6f}")
print(f"  alpha = {popt_int[1]:.4f}")

# Compare
print(f"\nDifference in alpha: {abs(popt_simp[1] - popt_int[1]):.4f}")
print(f"Ratio of C values: {popt_simp[0]/popt_int[0]:.4f}")

# Also check: what IS the integral value?
print("\n" + "=" * 72)
print("INTEGRAL VALUES (for reference)")
print("=" * 72)
for x in [1e6, 1e8, 1e10, 1e12, 1e14]:
    simp = hl_simplified(x)
    integ = hl_integral(x)
    diff_pct = (simp - integ) / integ * 100
    print(f"x = {x:.0e}: simplified = {simp:.1f}, integral = {integ:.1f}, diff = {diff_pct:.2f}%")
