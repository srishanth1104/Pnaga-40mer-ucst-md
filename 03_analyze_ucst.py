import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TEMPS = [30, 40, 50, 60, 70, 80]
OUT_DIR = "analysis"
os.makedirs(OUT_DIR, exist_ok=True)


def read_xvg(path, col=1):
    t, y = [], []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("#", "@")) or not line.strip():
                continue
            parts = line.split()
            t.append(float(parts[0]))
            y.append(float(parts[col]))
    return np.array(t), np.array(y)


rows = []
rg_series = {}
qc_series = {}

for T in TEMPS:
    rdir = os.path.join("runs", f"T{T}C")
    rg_path = os.path.join(rdir, "gyrate.xvg")
    if not os.path.exists(rg_path):
        print(f"Skipping {T} C (gyrate.xvg not found)")
        continue

    t_rg, rg = read_xvg(rg_path, col=1)
    rg_series[T] = (t_rg / 1000.0, rg)  # convert ps -> ns
    half = len(rg) // 2
    row = {"Temp_C": T, "Rg_mean_nm": rg[half:].mean(), "Rg_std_nm": rg[half:].std()}

    for key, fname in [
        ("SASA", "sasa.xvg"),
        ("IntraHB", "hbond_intra.xvg"),
        ("WaterHB", "hbond_water.xvg"),
    ]:
        fpath = os.path.join(rdir, fname)
        if os.path.exists(fpath):
            _, vals = read_xvg(fpath, col=1)
            row[f"{key}_mean"] = vals[len(vals) // 2 :].mean()
            row[f"{key}_std"] = vals[len(vals) // 2 :].std()

    rows.append(row)

    # Load NPT equilibration QC if present
    qc_entry = {}
    for qkey, qfile in [("T", "npt_temp.xvg"), ("D", "npt_dens.xvg")]:
        qpath = os.path.join(rdir, qfile)
        if os.path.exists(qpath):
            qc_entry[qkey] = read_xvg(qpath, col=1)
    if qc_entry:
        qc_series[T] = qc_entry

if not rows:
    raise SystemExit("No completed gyrate.xvg files found under runs/T*C/.")

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT_DIR, "ucst_summary.csv"), index=False)
print(df.to_string(index=False))

# 1. Main 4-panel UCST transition figure
fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=200)
axes = axes.ravel()

axes[0].errorbar(df["Temp_C"], df["Rg_mean_nm"], yerr=df["Rg_std_nm"], fmt="o-", color="#0F6E6E", capsize=4)
axes[0].set_ylabel("Radius of Gyration, Rg (nm)")
axes[0].set_title("Chain Expansion (UCST Transition)")

if "SASA_mean" in df.columns:
    axes[1].errorbar(df["Temp_C"], df["SASA_mean"], yerr=df["SASA_std"], fmt="s-", color="#1D4ED8", capsize=4)
    axes[1].set_ylabel("SASA (nm²)")
    axes[1].set_title("Solvent Accessible Surface Area")

if "IntraHB_mean" in df.columns:
    axes[2].errorbar(df["Temp_C"], df["IntraHB_mean"], yerr=df["IntraHB_std"], fmt="^-", color="#C2591A", capsize=4)
    axes[2].set_ylabel("Intra-Chain H-Bonds")
    axes[2].set_title("Polymer-Polymer H-Bonds")

if "WaterHB_mean" in df.columns:
    axes[3].errorbar(df["Temp_C"], df["WaterHB_mean"], yerr=df["WaterHB_std"], fmt="d-", color="#0369A1", capsize=4)
    axes[3].set_ylabel("Polymer-Water H-Bonds")
    axes[3].set_title("Polymer Hydration H-Bonds")

for ax in axes:
    ax.set_xlabel("Temperature (°C)")
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pnaga40_ucst_curve.png"))
plt.close()

# 2. Rg(t) time-series convergence plot
plt.figure(figsize=(9, 4.5), dpi=200)
for T, (t_ns, rg) in sorted(rg_series.items()):
    plt.plot(t_ns, rg, lw=1.2, alpha=0.85, label=f"{T} °C")
plt.xlabel("Simulation Time (ns)")
plt.ylabel("Radius of Gyration, Rg (nm)")
plt.title("40-mer PNAGA Rg Convergence Across Temperatures")
plt.legend(ncol=3, frameon=False)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "rg_timeseries.png"))
plt.close()

# 3. Equilibration QC plot (NPT temperature and density)
if qc_series:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), dpi=200)
    for T, qd in sorted(qc_series.items()):
        if "T" in qd:
            ax1.plot(qd["T"][0], qd["T"][1], lw=0.9, alpha=0.8, label=f"{T} °C")
        if "D" in qd:
            ax2.plot(qd["D"][0], qd["D"][1], lw=0.9, alpha=0.8, label=f"{T} °C")
    ax1.set_xlabel("NPT Time (ps)")
    ax1.set_ylabel("Temperature (K)")
    ax1.set_title("NPT Temperature Stability")
    ax2.set_xlabel("NPT Time (ps)")
    ax2.set_ylabel("Density (kg/m³)")
    ax2.set_title("NPT Density Equilibration")
    for ax in (ax1, ax2):
        ax.legend(ncol=2, fontsize=8, frameon=False)
        ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "equilibration_qc.png"))
    plt.close()

print("Saved all analysis outputs in analysis/")
