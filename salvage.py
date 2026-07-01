"""
CAN WE SALVAGE THIS?
Let's find what DOES follow a power law with the full integral.
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
    result, _ = quad(lambda t: 1/math.log(t), 2, x)
    return result

def hl_simplified(x):
    return 2 * C2 * x / math.log(x)**2

def hl_integral(x):
    return 2 * C2 * (li(x) - x/math.log(x))

def power_law(x, C, alpha):
    return C * x**alpha

xs = np.array(sorted(nicely_data.keys()))
pi2s = np.array([nicely_data[int(x)] for x in xs])

# ================================================================
# CANDIDATE 1: The correction term Delta(x) = HL_integral - HL_simplified
# This is what you need to ADD to the simplified to get the integral
# ================================================================
print("=" * 72)
print("CANDIDATE 1: Correction term Delta(x) = HL_integral - HL_simplified")
print("=" * 72)

Delta = np.array([hl_integral(x) - hl_simplified(x) for x in xs])
popt_d, _ = curve_fit(power_law, xs, np.abs(Delta), p0=[1, 0.5], maxfev=10000)
pred_d = power_law(xs, *popt_d)
r2_d = 1 - np.sum((np.abs(Delta) - pred_d)**2)/np.sum((np.abs(Delta) - np.mean(np.abs(Delta)))**2)
print(f"Delta(x) = {popt_d[0]:.4e} * x^{popt_d[1]:.4f}")
print(f"R^2 = {r2_d:.6f}")
print(f"alpha = {popt_d[1]:.4f}")

# ================================================================
# CANDIDATE 2: Relative error = Delta(x) / HL_integral(x)
# ================================================================
print(f"\n{'=' * 72}")
print("CANDIDATE 2: Relative error = Delta / HL_integral")
print(f"{'=' * 72}")

rel_err = np.array([(hl_integral(xs[i]) - hl_simplified(xs[i])) / hl_integral(xs[i]) for i in range(len(xs))])
print("Relative error values:")
for i in range(0, len(xs), 5):
    print(f"  x = {xs[i]:.0e}: {rel_err[i]*100:.2f}%")

# ================================================================
# CANDIDATE 3: R_int / x^{1/2} (normalized residual)
# ================================================================
print(f"\n{'=' * 72}")
print("CANDIDATE 3: R_int / x^{1/2}")
print(f"{'=' * 72}")

R_int = np.array([pi2s[i] - hl_integral(xs[i]) for i in range(len(xs))])
normalized = R_int / np.sqrt(xs)
popt_n, _ = curve_fit(power_law, xs, np.abs(normalized), p0=[1, 0.3], maxfev=10000)
pred_n = power_law(xs, *popt_n)
r2_n = 1 - np.sum((np.abs(normalized) - pred_n)**2)/np.sum((np.abs(normalized) - np.mean(np.abs(normalized)))**2)
print(f"|R_int|/sqrt(x) = {popt_n[0]:.4e} * x^{popt_n[1]:.4f}")
print(f"R^2 = {r2_n:.6f}")

# ================================================================
# CANDIDATE 4: Li(x) - x/log(x) - x/(log x)^2
# The "second correction" in the asymptotic expansion of Li(x)
# ================================================================
print(f"\n{'=' * 72}")
print("CANDIDATE 4: Li(x) - x/log(x) - x/(log x)^2")
print(f"{'=' * 72}")

def second_correction(x):
    return li(x) - x/math.log(x) - x/math.log(x)**2

corrections = np.array([second_correction(x) for x in xs])
popt_c, _ = curve_fit(power_law, xs, corrections, p0=[1, 0.5], maxfev=10000)
pred_c = power_law(xs, *popt_c)
r2_c = 1 - np.sum((corrections - pred_c)**2)/np.sum((corrections - np.mean(corrections))**2)
print(f"Li(x) - x/log(x) - x/(log x)^2 = {popt_c[0]:.4e} * x^{popt_c[1]:.4f}")
print(f"R^2 = {r2_c:.6f}")
print(f"alpha = {popt_c[1]:.4f}")

# ================================================================
# CANDIDATE 5: The "exact" HL prediction using asymptotic expansion
# HL_exact = 2C2 * (Li(x) - x/log(x) - x/(log x)^2 - 2x/(log x)^3)
# ================================================================
print(f"\n{'=' * 72}")
print("CANDIDATE 5: Higher-order HL correction")
print(f"{'=' * 72}")

def hl_exact(x):
    logx = math.log(x)
    return 2 * C2 * (li(x) - x/logx - x/logx**2 - 2*x/logx**3)

R_exact = np.array([pi2s[i] - hl_exact(xs[i]) for i in range(len(xs))])
popt_e, _ = curve_fit(power_law, xs, np.abs(R_exact), p0=[1, 0.5], maxfev=10000)
pred_e = power_law(xs, *popt_e)
r2_e = 1 - np.sum((np.abs(R_exact) - pred_e)**2)/np.sum((np.abs(R_exact) - np.mean(np.abs(R_exact)))**2)
print(f"|R_exact| = {popt_e[0]:.4e} * x^{popt_e[1]:.4f}")
print(f"R^2 = {r2_e:.6f}")

# ================================================================
# CANDIDATE 6: pi_2(x) - 2C2 * Li(x)  (subtracting Li directly)
# ================================================================
print(f"\n{'=' * 72}")
print("CANDIDATE 6: pi_2(x) - 2C2 * Li(x)")
print(f"{'=' * 72}")

R_li = np.array([pi2s[i] - 2*C2*li(xs[i]) for i in range(len(xs))])
popt_l, _ = curve_fit(power_law, xs, np.abs(R_li), p0=[1, 0.5], maxfev=10000)
pred_l = power_law(xs, *popt_l)
r2_l = 1 - np.sum((np.abs(R_li) - pred_l)**2)/np.sum((np.abs(R_li) - np.mean(np.abs(R_li)))**2)
print(f"|R_Li| = {popt_l[0]:.4e} * x^{popt_l[1]:.4f}")
print(f"R^2 = {r2_l:.6f}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'=' * 72}")
print("SUMMARY")
print(f"{'=' * 72}")
print(f"{'Candidate':<45} {'alpha':>8} {'R^2':>8}")
print("-" * 63)
print(f"{'Delta(x) = HL_int - HL_simp':<45} {popt_d[1]:8.4f} {r2_d:8.4f}")
print(f"{'Li(x) - x/log(x) - x/(log x)^2':<45} {popt_c[1]:8.4f} {r2_c:8.4f}")
print(f"{'|R_exact| (3-term correction)':<45} {popt_e[1]:8.4f} {r2_e:8.4f}")
print(f"{'|R_Li| = |pi_2 - 2C2*Li(x)|':<45} {popt_l[1]:8.4f} {r2_l:8.4f}")
