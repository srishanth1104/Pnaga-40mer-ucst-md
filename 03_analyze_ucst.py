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
for T in TEMPS:
    rdir = os.path.join("runs", f"T{T}C")
    rg_path = os.path.join(rdir, "gyrate.xvg")
    if not os.path.exists(rg_path):
        print(f"Skipping {T} C (gyrate.xvg not found)")
        continue
    _, rg = read_xvg(rg_path, col=1)
    half = len(rg) // 2
    row = {"Temp_C": T, "Rg_mean_nm": rg[half:].mean(), "Rg_std_nm": rg[half:].std()}

    sasa_path = os.path.join(rdir, "sasa.xvg")
    if os.path.exists(sasa_path):
        _, sasa = read_xvg(sasa_path, col=1)
        row["SASA_mean_nm2"] = sasa[len(sasa) // 2 :].mean()
        row["SASA_std_nm2"] = sasa[len(sasa) // 2 :].std()

    hb_path = os.path.join(rdir, "hbond_intra.xvg")
    if os.path.exists(hb_path):
        _, hb = read_xvg(hb_path, col=1)
        row["IntraHB_mean"] = hb[len(hb) // 2 :].mean()
        row["IntraHB_std"] = hb[len(hb) // 2 :].std()

    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT_DIR, "ucst_summary.csv"), index=False)
print(df.to_string(index=False))

fig, axes = plt.subplots(1, 3, figsize=(13, 4), dpi=200)
axes[0].errorbar(df["Temp_C"], df["Rg_mean_nm"], yerr=df["Rg_std_nm"], fmt="o-", color="#0F6E6E", capsize=4)
axes[0].set_xlabel("Temperature (°C)")
axes[0].set_ylabel("Radius of Gyration, Rg (nm)")
axes[0].set_title("Chain Expansion (UCST)")

if "SASA_mean_nm2" in df.columns:
    axes[1].errorbar(df["Temp_C"], df["SASA_mean_nm2"], yerr=df["SASA_std_nm2"], fmt="s-", color="#1D4ED8", capsize=4)
    axes[1].set_xlabel("Temperature (°C)")
    axes[1].set_ylabel("SASA (nm²)")
    axes[1].set_title("Solvent Exposure")

if "IntraHB_mean" in df.columns:
    axes[2].errorbar(df["Temp_C"], df["IntraHB_mean"], yerr=df["IntraHB_std"], fmt="^-", color="#C2591A", capsize=4)
    axes[2].set_xlabel("Temperature (°C)")
    axes[2].set_ylabel("Intramolecular H-Bonds")
    axes[2].set_title("Polymer-Polymer H-Bonds")

for ax in axes:
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pnaga40_ucst_curve.png"))
print("Saved analysis/pnaga40_ucst_curve.png and analysis/ucst_summary.csv")
