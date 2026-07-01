"""
Extended power law analysis with Nicely's full dataset.
Tests theoretical connections, bootstrap CI, and extended data.
"""
import math
import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy import stats

C2 = 0.6601618914

# Nicely's verified data: pi_2(x) for x = 10^6 to 10^14
# Sources: Nicely (1995, 1999), verified tables
nicely_data = {
    10**6:  8169,
    10**7:  58980,
    10**8:  440312,
    10**9:  3424506,
    10**10: 27412679,
    10**11: 224376048,
    # Extended data from Nicely's tables
    2e11:   424084653,
    3e11:   615885700,
    4e11:   802817718,
    5e11:   986222314,
    6e11:   1166916933,
    7e11:   1345394380,
    8e11:   1521998439,
    9e11:   1696987738,
    10**12: 1870585220,
    2e12:   3552770943,
    3e12:   5173760785,
    4e12:   6756832076,
    5e12:   8312493003,
    6e12:   9846842484,
    7e12:   11363874338,
    8e12:   12866256870,
    9e12:   14356002120,
    10**13: 15834664872,
    2e13:   30198862775,
    3e13:   44078684643,
    4e13:   57657248284,
    5e13:   71018282471,
    6e13:   84209699420,
    7e13:   97262712867,
    8e13:   110198743491,
    9e13:   123033833767,
    10**14: 135780321665,
}

def hl(x):
    return 2 * C2 * x / math.log(x)**2

def residual(x):
    return nicely_data[int(x)] - hl(x)

def power_law(x, C, alpha):
    return C * x**alpha

# Sort data
xs = sorted(nicely_data.keys())
R = [residual(x) for x in xs]

print("=" * 70)
print("EXTENDED POWER LAW ANALYSIS — Nicely's Full Dataset")
print("=" * 70)
print(f"\nData points: {len(xs)} (from 10^6 to 10^14)")
print(f"Range: {xs[0]:.0e} to {xs[-1]:.0e}")

# === 1. POWER LAW FIT ON EXTENDED DATA ===
print("\n" + "=" * 70)
print("1. POWER LAW FIT ON EXTENDED DATA")
print("=" * 70)

log_x = np.log10(xs)
log_R = np.log10(R)

# Fit all points
coeffs = np.polyfit(log_x, log_R, 1)
alpha_full = coeffs[0]
C_full = 10**coeffs[1]
print(f"\nFull fit ({len(xs)} points): R(x) = {C_full:.4e} * x^{alpha_full:.4f}")

# Fit residuals
preds = [power_law(x, C_full, alpha_full) for x in xs]
ss_res = sum((p - r)**2 for p, r in zip(preds, R))
ss_tot = sum((r - np.mean(R))**2 for r in R)
r2 = 1 - ss_res / ss_tot
rmse = math.sqrt(ss_res)
print(f"R^2 = {r2:.10f}")
print(f"RMSE = {rmse:.1f}")
print(f"Mean relative error = {np.mean([abs(p-r)/abs(r)*100 for p, r in zip(preds, R)]):.1f}%")

# === 2. LOCAL EXPONENT ANALYSIS ===
print("\n" + "=" * 70)
print("2. LOCAL EXPONENT DRIFT")
print("=" * 70)

print(f"\n{'x_i':>12s}  {'x_{i+1}':>12s}  {'alpha_local':>12s}  {'C_local':>12s}")
local_alphas = []
for i in range(len(xs) - 1):
    x1, x2 = xs[i], xs[i+1]
    r1, r2 = R[i], R[i+1]
    alpha_loc = math.log(r2/r1) / math.log(x2/x1)
    C_loc = r1 / x1**alpha_loc
    local_alphas.append(alpha_loc)
    print(f"  {x1:>10.0e}  {x2:>10.0e}  {alpha_loc:>12.4f}  {C_loc:>12.4e}")

print(f"\nAlpha drift: {min(local_alphas):.4f} to {max(local_alphas):.4f}")
print(f"Mean alpha: {np.mean(local_alphas):.4f} ± {np.std(local_alphas):.4f}")

# === 3. BOOTSTRAP CONFIDENCE INTERVALS ===
print("\n" + "=" * 70)
print("3. BOOTSTRAP CONFIDENCE INTERVALS ON ALPHA")
print("=" * 70)

n_bootstrap = 10000
bootstrap_alphas = []
bootstrap_Cs = []

for _ in range(n_bootstrap):
    indices = np.random.choice(len(xs), size=len(xs), replace=True)
    bx = [xs[i] for i in indices]
    bR = [R[i] for i in indices]
    log_bx = np.log10(bx)
    log_bR = np.log10(bR)
    bc = np.polyfit(log_bx, log_bR, 1)
    bootstrap_alphas.append(bc[0])
    bootstrap_Cs.append(10**bc[1])

alpha_ci = np.percentile(bootstrap_alphas, [2.5, 50, 97.5])
C_ci = np.percentile(bootstrap_Cs, [2.5, 50, 97.5])
print(f"\nAlpha: {alpha_ci[1]:.4f} (95% CI: [{alpha_ci[0]:.4f}, {alpha_ci[2]:.4f}])")
print(f"C:     {C_ci[1]:.4e} (95% CI: [{C_ci[0]:.4e}, {C_ci[2]:.4e}])")

# === 4. CONNECTION TO KNOWN CONSTANTS ===
print("\n" + "=" * 70)
print("4. CONNECTION TO KNOWN MATHEMATICAL CONSTANTS")
print("=" * 70)

candidates = {
    "5/6": 5/6,
    "euler-mascheroni": 0.5772156649,
    "Catalan (G)": 0.9159655942,
    "1 - 1/e": 1 - 1/math.e,
    "ln(2)": math.log(2),
    "pi/4": math.pi/4,
    "sqrt(2)-1": math.sqrt(2)-1,
    "1/sqrt(2*pi)": 1/math.sqrt(2*math.pi),
    "zeta(3)/zeta(2)": 1.2020569/1.6449341,
    "2/pi": 2/math.pi,
    "1/sqrt(3)": 1/math.sqrt(3),
    "ln(pi)": math.log(math.pi),
}

print(f"\nFitted alpha: {alpha_full:.6f}")
print(f"\n{'Constant':<25s}  {'Value':>12s}  {'dist':>12s}  {'sigma away':>10s}")
print("-" * 65)

best_name = None
best_dist = float('inf')
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - alpha_full)):
    dist = abs(alpha_full - val)
    sigma_dist = dist / np.std(bootstrap_alphas)
    marker = " <-- CLOSEST" if dist < best_dist else ""
    if dist < best_dist:
        best_dist = dist
        best_name = name
    print(f"  {name:<25s}  {val:>12.6f}  {dist:>12.6f}  {sigma_dist:>8.1f}sigma{marker}")

print(f"\nClosest: {best_name} (distance = {best_dist:.6f}, {best_dist/np.std(bootstrap_alphas):.1f}sigma)")

# === 5. LOG-LOG CURVATURE TEST ===
print("\n" + "=" * 70)
print("5. IS ALPHA TRULY CONSTANT? (Curvature test)")
print("=" * 70)

# Fit quadratic in log-log space: log(R) = a + b*log(x) + c*log(x)^2
log_x_arr = np.array(log_x)
log_R_arr = np.array(log_R)
coeffs_quad = np.polyfit(log_x_arr, log_R_arr, 2)
c_quad = coeffs_quad[0]

print(f"\nQuadratic fit: log(R) = {coeffs_quad[2]:.4f} + {coeffs_quad[1]:.4f}*log(x) + {coeffs_quad[0]:.4f}*log(x)^2")
print(f"Quadratic term: c = {c_quad:.6f}")

if abs(c_quad) < 0.001:
    print("=> Alpha is essentially constant (no curvature)")
elif c_quad > 0:
    print("=> Alpha is slowly INCREASING (positive curvature)")
else:
    print("=> Alpha is slowly DECREASING (negative curvature)")

# Compute effective alpha at each scale
eff_alphas = []
for i in range(len(xs)):
    log_x_arr_i = log_x_arr.copy()
    log_R_arr_i = log_R_arr.copy()
    # Local quadratic fit
    if i < 2:
        fit_range = slice(0, min(5, len(xs)))
    elif i > len(xs) - 3:
        fit_range = slice(max(0, len(xs)-5), len(xs))
    else:
        fit_range = slice(i-2, i+3)
    local_coeffs = np.polyfit(log_x_arr_i[fit_range], log_R_arr_i[fit_range], 2)
    eff_alpha = local_coeffs[1] + 2 * local_coeffs[0] * log_x_arr[i]
    eff_alphas.append(eff_alpha)

print(f"\nEffective alpha at each scale:")
for x, ea in zip(xs, eff_alphas):
    print(f"  x={x:>10.0e}  alpha_eff={ea:.4f}")

# === 6. PREDICTION FOR 10^15 ===
print("\n" + "=" * 70)
print("6. PREDICTIONS FOR LARGER SCALES")
print("=" * 70)

for x_pred in [10**15, 10**16]:
    R_pred = power_law(x_pred, C_full, alpha_full)
    hl_val = hl(x_pred)
    pi2_pred = hl_val + R_pred
    print(f"\n  x = {x_pred:.0e}:")
    print(f"    R(x) = {R_pred:.2e}")
    print(f"    pi_2(x) = {pi2_pred:.0e}")
    
    # Extrapolation uncertainty
    R_low = power_law(x_pred, C_ci[0], alpha_ci[0])
    R_high = power_law(x_pred, C_ci[2], alpha_ci[2])
    pi2_low = hl(x_pred) + R_low
    pi2_high = hl(x_pred) + R_high
    print(f"    95% CI: [{pi2_low:.0e}, {pi2_high:.0e}]")

# === 7. GOODNESS OF FIT ===
print("\n" + "=" * 70)
print("7. GOODNESS OF FIT STATISTICS")
print("=" * 70)

# Chi-squared test (assuming Poisson errors on pi_2)
chi2 = 0
for x, r, p in zip(xs, R, preds):
    pi2_val = nicely_data[int(x)]
    # Variance of pi_2 is approximately pi_2 (Poisson)
    var_r = pi2_val  # since R = pi_2 - HL, and HL is deterministic
    chi2 += (r - p)**2 / var_r

dof = len(xs) - 2  # 2 parameters
p_value = 1 - stats.chi2.cdf(chi2, dof)
print(f"\nChi-squared: {chi2:.2f}")
print(f"Degrees of freedom: {dof}")
print(f"P-value: {p_value:.6f}")
print(f"Reduced chi-squared: {chi2/dof:.4f}")

# === 8. SUMMARY ===
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
Extended analysis with {len(xs)} data points (10^6 to 10^14):

1. Power law R(x) = C * x^alpha holds with:
   - alpha = {alpha_full:.4f} +/- {np.std(bootstrap_alphas):.4f} (95% CI: [{alpha_ci[0]:.4f}, {alpha_ci[2]:.4f}])
   - C = {C_full:.4e}
   - R^2 = {r2:.8f}

2. Alpha drift: {min(local_alphas):.4f} to {max(local_alphas):.4f}
   {'Alpha is increasing - possible log correction needed' if max(local_alphas) - min(local_alphas) > 0.05 else 'Alpha is essentially constant'}

3. Closest known constant: {best_name} (distance: {best_dist:.6f}, {best_dist/np.std(bootstrap_alphas):.1f}sigma)

4. Quadratic curvature term: c = {c_quad:.6f}
   {'No significant curvature' if abs(c_quad) < 0.001 else 'Significant curvature detected'}

5. Chi-squared p-value: {p_value:.4f}
   {'Good fit' if p_value > 0.05 else 'Poor fit - power law may not be exact'}
""")
