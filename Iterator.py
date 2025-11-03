# ================================================================
# chord_span_taper_sweep_with_pitch_plots.py
# ================================================================
# Sweeps chord (root), span, and taper ratio for a trapezoidal wing
# Uses NACA 4415 polars and computes aerodynamic + pitch-up metrics.
# ================================================================

import os, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------
# INPUT POLAR FILES
# -----------------------
POLAR_500K = "xf-naca4415-il-500000.csv"
POLAR_1M   = "xf-naca4415-il-1000000.csv"

# -----------------------
# AIRCRAFT & ENVIRONMENT BASED ON THE Aerosonde Mk 4.7J and Boeing Insitu ScanEagle
# -----------------------
W = 220          # N
P_engine = 5000   # W (shaft power)
eta = 0.9
m = 1.0
alt_cruise = 2743.2  # m
CD0 = 0.0248

# -----------------------
# ATMOSPHERE FUNCTIONS
# -----------------------
def isa_temp(alt):
    T0 = 288.16; lapse = -0.0065
    return T0 + lapse * min(alt, 11000)

def isa_density(alt):
    T0 = 288.16; p0 = 101325.0; g = 9.80665; R = 287.05; lapse = -0.0065
    if alt <= 11000:
        T = T0 + lapse * alt
        p = p0 * (T / T0) ** (-g / (lapse * R))
        return p / (R * T)
    T11 = T0 + lapse * 11000
    p11 = p0 * (T11 / T0) ** (-g / (lapse * R))
    H = R * T11 / g
    rho11 = p11 / (R * T11)
    return rho11 * math.exp(-(alt - 11000) / H)

def mu_sutherland(T):
    mu_ref = 1.7894e-5
    T_ref = 288.15
    S = 110.4
    return mu_ref * (T / T_ref)**1.5 * (T_ref + S) / (T + S)

rho0 = isa_density(0)
rho = isa_density(alt_cruise)
T = isa_temp(alt_cruise)
mu = mu_sutherland(T)

# -----------------------
# LOAD POLARS
# -----------------------
def load_polar(path):
    header_idx = None
    with open(path, 'r', errors='ignore') as fh:
        for i, line in enumerate(fh):
            if line.strip().startswith("Alpha"):
                header_idx = i; break
    if header_idx is None:
        raise ValueError(f"Header not found in {path}")
    df = pd.read_csv(path, skiprows=header_idx)
    df.columns = [c.strip() for c in df.columns]
    return df

def find_clcol(df):
    for c in df.columns:
        if c.lower().startswith('cl'):
            return c
    return df.columns[1]

df500 = load_polar(POLAR_500K)
df1m  = load_polar(POLAR_1M)
clcol500, clcol1m = find_clcol(df500), find_clcol(df1m)
CL2D_500k = float(df500[clcol500].max())
CL2D_1M   = float(df1m[clcol1m].max())
Re_low, Re_high = 5e5, 1e6

# -----------------------
# SWEEP PARAMETERS
# -----------------------
root_chords = np.linspace(0, 1, 25)  # m
span_vals   = np.linspace(0, 3.5, 25)   # m
taper_vals  = np.linspace(0, 1, 25)    # c_t / c_r

# -----------------------
# HELPER FUNCTIONS
# -----------------------
def span_efficiency_from_taper(lam):
    # approximate curve peaking near λ=0.6
    return 0.9 * (1 - 0.08 * (lam - 0.4)**2)

def mean_aero_chord(c_r, lam):
    return (2/3)*c_r*((1+lam+lam**2)/(1+lam))

def K_from_AR(AR, e): return 1.0 / (math.pi * e * AR)

results = []
for c_r in root_chords:
    for b in span_vals:
        for lam in taper_vals:
            c_t = lam * c_r
            S = 0.5 * b * (c_r + c_t)
            AR = b**2 / S
            e = span_efficiency_from_taper(lam)
            K = K_from_AR(AR, e)
            c_mac = mean_aero_chord(c_r, lam)
            inner = K / (3.0 * CD0)
            if inner <= 0: continue
            Vend = math.sqrt((2.0 * W) / (rho * S)) * inner**0.25
            if Vend <= 0: continue
            Re = rho * Vend * c_mac / mu

            # interpolate CL2D
            logR = math.log10(Re)
            log_low, log_high = math.log10(Re_low), math.log10(Re_high)
            if Re <= Re_low:
                CL2 = CL2D_500k
            elif Re >= Re_high:
                CL2 = CL2D_1M
            else:
                t = (logR - log_low) / (log_high - log_low)
                CL2 = CL2D_500k + t*(CL2D_1M - CL2D_500k)

            CL3 = CL2 / (1 + CL2 / (math.pi * e * AR))
            CL_req = W / (0.5 * rho * Vend**2 * S)
            if CL_req > CL3: continue
            Vstall = math.sqrt((2.0 * W) / (rho * S * CL3))
            endurance_ratio = Vend / Vstall
            if endurance_ratio < 0: continue

            Dp = 0.5 * rho * Vend**2 * S * CD0
            Di = 2 * K * W**2 / (rho * Vend**2 * S)
            P_req = (Dp + Di) * Vend
            P_avail = P_engine * eta * (rho / rho0)**m
            P_excess = P_avail - P_req
            score = -P_req/1000 + 0.4*endurance_ratio

            results.append({
                "c_root": c_r, "c_tip": c_t, "taper": lam,
                "b": b, "S": S, "AR": AR, "e": e, "Re": Re,
                "CL3D": CL3, "CL_req": CL_req, "Vend": Vend, "Vstall": Vstall,
                "EnduranceRatio": endurance_ratio,
                "P_req_kW": P_req/1000, "P_excess_kW": P_excess/1000,
                "Score": score
            })

df = pd.DataFrame(results)
df.sort_values("Score", ascending=False, inplace=True)
df.reset_index(drop=True, inplace=True)

# -----------------------
# PITCH-UP METRICS
# -----------------------
def estimate_a2D(df):
    a_col = df.columns[0]; cl_col = find_clcol(df)
    df_lin = df[(df[a_col] >= -2) & (df[a_col] <= 8)]
    if len(df_lin) < 3: return 2 * math.pi
    coeffs = np.polyfit(df_lin[a_col], df_lin[cl_col], 1)
    a_per_deg = coeffs[0]
    return a_per_deg * (180 / np.pi)

a2D_500k, a2D_1M = estimate_a2D(df500), estimate_a2D(df1m)

def compute_pitch(row):
    c_mac, AR, CL3, CLreq, lam = mean_aero_chord(row["c_root"], row["taper"]), row["AR"], row["CL3D"], row["CL_req"], row["taper"]
    e = span_efficiency_from_taper(lam)
    Re_here = rho * row["Vend"] * c_mac / mu
    logR = math.log10(Re_here)
    log_low, log_high = math.log10(Re_low), math.log10(Re_high)
    if Re_here <= Re_low: a2D = a2D_500k
    elif Re_here >= Re_high: a2D = a2D_1M
    else:
        t = (logR - log_low) / (log_high - log_low)
        a2D = a2D_500k + t*(a2D_1M - a2D_500k)
    a3D = a2D / (1 + a2D / (math.pi * e * AR))
    dCL = max(CL3 - CLreq, 0)
    dAlpha_deg = math.degrees(dCL / a3D)
    n_max = CL3 / CLreq if CLreq > 0 else np.nan
    return pd.Series({"a3D_per_rad": a3D, "delta_alpha_deg": dAlpha_deg, "n_max": n_max})

df = pd.concat([df, df.apply(compute_pitch, axis=1)], axis=1)

print("\nTop 10 designs:")
print(df.head(10).to_string(index=False))

# -----------------------
# PLOTTING
# -----------------------
def plot_group(df, xvar, varlabel):
    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(f"Trends vs {varlabel}")
    axs = axs.flatten()

    def add(ax, yvar, cmap, label):
        sc = ax.scatter(df[xvar], df[yvar], c=df["Score"], cmap=cmap)
        ax.set_xlabel(varlabel); ax.set_ylabel(label); ax.grid(True)
        cbar = plt.colorbar(sc, ax=ax); cbar.set_label("Score")

    add(axs[0], "P_req_kW", "jet_r", "Power Req (kW)")
    add(axs[1], "EnduranceRatio", "jet_r", "Vend / Vstall")
    add(axs[2], "Vstall", "jet_r", "Vstall (m/s)")
    add(axs[3], "Score", "jet_r", "Score")
    add(axs[4], "delta_alpha_deg", "jet_r", "Δα before stall (deg)")
    add(axs[5], "n_max", "jet_r", "n_max (Load Factor)")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig

figs = [
    plot_group(df, "AR", "Aspect Ratio"),
    plot_group(df, "c_root", "Root Chord (m)"),
    plot_group(df, "b", "Wingspan (m)"),
    plot_group(df, "S", "Planform Area (m²)"),
    plot_group(df, "taper", "Taper Ratio (λ)")
]
plt.show(block=True)
