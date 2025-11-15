# ================================================================
# chord_span_sweep_with_pitch_plots_v2.py
# ================================================================
# Sweeps chord & span for rectangular wing using NACA 4415 polars
# Adds aerodynamic + pitch-up metrics, and produces simultaneous
# plots for AR, chord, span, and planform area.
# ================================================================

import os, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------
# INPUT POLAR FILES
# -----------------------
#POLAR_500K = "xf-naca4415-il-500000.csv"
#POLAR_1M   = "xf-naca4415-il-1000000.csv"
POLAR_500K = "xf-sg6043-il-500000.csv"
POLAR_1M   = "xf-sg6043-il-1000000.csv"

# -----------------------
# AIRCRAFT & ENVIRONMENT
# -----------------------
W = 294.3          # N
P_engine = 3000   # W (shaft power)
eta = 0.9         # prop efficiency
m = 1.0           # power altitude exponent
alt_cruise = 2026.92  # m
CD0 = 0.0248
e = .7           # span efficiency factor

# -----------------------
# STANDARD ATMOSPHERE
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
    else:
        T11 = T0 + lapse * 11000
        p11 = p0 * (T11 / T0) ** (-g / (lapse * R))
        H = R * T11 / g
        rho11 = p11 / (R * T11)
        return rho11 * math.exp(-(alt - 11000) / H)

def mu_sutherland(T):
    mu_ref = 1.7894e-5  # Pa·s at 288.15 K
    T_ref = 288.15
    S = 110.4
    return mu_ref * (T / T_ref) ** 1.5 * (T_ref + S) / (T + S)

rho0 = isa_density(0)
rho = isa_density(alt_cruise)
T = isa_temp(alt_cruise)
mu = mu_sutherland(T)

# -----------------------
# LOAD POLARS
# -----------------------
def load_polar(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
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
chord_vals = np.arange(0.25, .325, 0.025)
span_vals  = np.arange(0.10, 3.1, 0.1)

def K_from_AR(AR, e): return 1.0 / (math.pi * e * AR)

results = []
for c in chord_vals:
    for b in span_vals:
        S = b * c
        AR = b / c
        K = K_from_AR(AR, e)
        inner = K / (3.0 * CD0)
        Vend = 26.8224#math.sqrt((2.0 * W) / (rho * S)) * inner ** 0.25
        if Vend <= 0: continue
        Re = rho * Vend * c / mu
        # interpolate CL2D
        logR = math.log10(Re)
        log_low, log_high = math.log10(Re_low), math.log10(Re_high)
        if Re <= Re_low:
            CL2 = CL2D_500k
        elif Re >= Re_high:
            CL2 = CL2D_1M
        else:
            t = (logR - log_low) / (log_high - log_low)
            CL2 = CL2D_500k + t * (CL2D_1M - CL2D_500k)
        # finite wing correction
        # --- helper: estimate 2D lift curve slope (per rad) from polar ---
        def estimate_a2D_from_polar(df_polar, alpha_col=None, cl_col=None, deg_range=(-2, 8)):
            # pick columns if not provided
            if alpha_col is None:
                alpha_col = df_polar.columns[0]
            if cl_col is None:
                for c in df_polar.columns:
                    if c.lower().startswith("cl"):
                        cl_col = c
                        break
            df_lin = df_polar[(df_polar[alpha_col] >= deg_range[0]) & (df_polar[alpha_col] <= deg_range[1])]
            x_deg = pd.to_numeric(df_lin[alpha_col], errors='coerce').values
            y_cl = pd.to_numeric(df_lin[cl_col], errors='coerce').values
            valid = (~np.isnan(x_deg)) & (~np.isnan(y_cl))
            if valid.sum() < 3:
                return 2.0 * np.pi  # fallback (thin airfoil)
            # fit CL = a_per_deg * alpha_deg + b
            a_per_deg = np.polyfit(x_deg[valid], y_cl[valid], 1)[0]
            # convert to per-radian
            a_per_rad = a_per_deg * (180.0 / np.pi)
            return a_per_rad


        # estimate once (outside loop)
        a2D_500 = estimate_a2D_from_polar(df500)
        a2D_1M = estimate_a2D_from_polar(df1m)

        # inside sweep loop, after computing Re and CL2 (CL2 is CL_max_2D from the polar)
        # compute a2D appropriate for this Re (log interpolation)
        logR = math.log10(Re)
        log_low = math.log10(Re_low)
        log_high = math.log10(Re_high)
        if Re <= Re_low:
            a2D_here = a2D_500
        elif Re >= Re_high:
            a2D_here = a2D_1M
        else:
            t = (logR - log_low) / (log_high - log_low)
            a2D_here = a2D_500 + t * (a2D_1M - a2D_500)

        # finite-wing slope
        a3D = a2D_here / (1.0 + (a2D_here / (math.pi * e * AR)))
        if a3D <= 0:
            a3D = 2.0 * math.pi / (1.0 + 2.0 * math.pi / (math.pi * e * AR))  # fallback

        # convert CLmax_2D to CLmax_3D by scaling with slope ratio
        CL3 = CL2 * (a3D / a2D_here)

        CL_req = W / (0.5 * rho * Vend**2 * S)
        if CL_req > CL3: continue
        Vstall = math.sqrt((2 * W) / (rho * S * CL3))
        endurance_ratio = Vend / Vstall
        if endurance_ratio < 1: continue
        # power calcs
        Dp = 0.5 * rho * Vend**2 * S * CD0
        Di = 2 * K * W**2 / (rho * Vend**2 * S)
        P_req = (Dp + Di) * Vend
        P_avail = P_engine * eta * (rho / rho0)**m
        P_excess = P_avail - P_req
        score = -P_req / 1000 + 0.4 * endurance_ratio
        results.append({
            "chord_m": c, "span_m": b, "S_m2": S, "AR": AR, "Re": Re,
            "CL3D": CL3, "CL_req": CL_req,
            "Vend": Vend, "Vstall": Vstall, "EnduranceRatio": endurance_ratio,
            "P_req_kW": P_req / 1000, "P_excess_kW": P_excess / 1000,
            "Score": score
        })

df = pd.DataFrame(results)
if df.empty:
    raise RuntimeError("No feasible designs found – adjust sweep or constraints.")
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
    c, AR, Vend, CL3, CLreq = row["chord_m"], row["AR"], row["Vend"], row["CL3D"], row["CL_req"]
    Re_here = rho * Vend * c / mu
    logR = math.log10(Re_here)
    log_low, log_high = math.log10(Re_low), math.log10(Re_high)
    if Re_here <= Re_low: a2D = a2D_500k
    elif Re_here >= Re_high: a2D = a2D_1M
    else:
        t = (logR - log_low) / (log_high - log_low)
        a2D = a2D_500k + t * (a2D_1M - a2D_500k)
    a3D = a2D / (1 + a2D / (math.pi * e * AR))
    dCL = max(CL3 - CLreq, 0)
    dAlpha_deg = math.degrees(dCL / a3D)
    n_max = CL3 / CLreq if CLreq > 0 else np.nan
    return pd.Series({"a3D_per_rad": a3D, "delta_alpha_deg": dAlpha_deg, "n_max": n_max})

df = pd.concat([df, df.apply(compute_pitch, axis=1)], axis=1)

print("\nTop 10 designs:")
print(df.head(10).to_string(index=False))

# -----------------------
# PLOTTING FUNCTION
# -----------------------
def plot_group(df, xvar, varlabel):
    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(f"Trends vs {varlabel}")
    axs = axs.flatten()

    def add_scatter(ax, yvar, cmap, label):
        sc = ax.scatter(df[xvar], df[yvar], c=df["Score"], cmap=cmap)
        ax.set_xlabel(varlabel)
        ax.set_ylabel(label)
        ax.grid(True)
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label("Score")

    add_scatter(axs[0], "P_req_kW", 'viridis', "Power Req (kW)")
    add_scatter(axs[1], "EnduranceRatio", 'plasma', "Vend / Vstall")
    add_scatter(axs[2], "Vstall", 'cool', "Vstall (m/s)")
    add_scatter(axs[3], "Score", 'magma', "Score")
    add_scatter(axs[4], "delta_alpha_deg", 'cividis', "Δα before stall (deg)")
    add_scatter(axs[5], "n_max", 'inferno', "n_max (Load Factor)")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig

# -----------------------
# GENERATE PLOTS SIMULTANEOUSLY
# -----------------------
figs = [
    plot_group(df, "AR", "Aspect Ratio"),
    plot_group(df, "chord_m", "Chord Length (m)"),
    plot_group(df, "span_m", "Wingspan (m)"),
    plot_group(df, "S_m2", "Planform Area (m²)")
]
plt.show(block=True)
