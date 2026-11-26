# 40-mer PNAGA UCST Molecular Dynamics (GROMACS GPU)

Complete, self-contained GROMACS GPU workflow to study the Upper Critical Solution Temperature (UCST) coil-to-globule transition of a single 40-mer poly(*N*-acryloyl glycinamide) (PNAGA) chain in explicit TIP3P water across **30, 40, 50, 60, 70, and 80 °C** (100 ns production per temperature, GAFF2 + AM1-BCC charges).

## Repository Contents

| File | Purpose |
| :--- | :--- |
| `00_check_env.sh` | Pre-flight check for GPU, `gmx`/`gmx_mpi` CUDA support, and Python/AmberTools |
| `environment.yml` | Conda environment for polymer building & UCST plotting (`rdkit`, `ambertools`, `acpype`, `matplotlib`) |
| `01_build_pnaga40.py` | Builds the 40-mer PNAGA chain (GAFF2 + AM1-BCC) and solvates it in a 1.5 nm dodecahedron TIP3P water box |
| `02_run_ucst_gpu.sh` | Runs EM → 200 ps NVT → 1 ns NPT → 100 ns GPU production MD + extracts `Rg`, `SASA`, intra-chain & polymer–water H-bonds, and QC metrics |
| `submit_slurm.sh` | SLURM GPU job-array script (`--array=0-5` runs all 6 temperatures in parallel) |
| `submit_pbs.sh` | PBS/Torque GPU script (alternative for non-SLURM clusters) |
| `03_analyze_ucst.py` | Generates `analysis/ucst_summary.csv`, `pnaga40_ucst_curve.png`, `rg_timeseries.png`, and `equilibration_qc.png` |
| `04_package_results.sh` | Packs only the lightweight results (`analysis/`, `.xvg` tables, and polymer-only `.xtc`/`.gro`/`.tpr`) into `PNAGA40_UCST_RESULTS.tar.gz` |

## 5-Step Execution Guide

1. **Clone & create the builder/analysis environment (one-time):**
   ```bash
   git clone https://github.com/srishanth1104/Pnaga-40mer-ucst-md.git
   cd Pnaga-40mer-ucst-md
   conda env create -f environment.yml
   conda activate pnaga_md
   ```
   Load your cluster's GPU-enabled GROMACS module (e.g., `module load gromacs`).

2. **Verify tools & GPU (10 seconds):**
   ```bash
   bash 00_check_env.sh
   ```

3. **Build & solvate the master 40-mer (~1 minute):**
   ```bash
   python 01_build_pnaga40.py
   ```
   Creates `00_master_inputs/solv_master.gro` and `00_master_inputs/topol_master.top`.

4. **Run all 6 temperatures (30, 40, 50, 60, 70, 80 °C × 100 ns) on GPU:**
   - **SLURM (parallel across 6 GPUs):** Edit `#SBATCH --partition` and `module load` in `submit_slurm.sh`, then run:
     ```bash
     sbatch submit_slurm.sh
     ```
   - **PBS/Torque:** Edit queue/module lines in `submit_pbs.sh`, then run `qsub submit_pbs.sh`.
   - **Direct single-GPU run (all 6 sequentially):**
     ```bash
     bash 02_run_ucst_gpu.sh
     ```
   - **Single temperature (e.g., 30 °C):**
     ```bash
     bash 02_run_ucst_gpu.sh 30
     ```
   *Auto-resume is enabled (`-cpi prod.cpt -cpt 15`): if a job hits a walltime limit, simply re-run the exact same command to resume cleanly.*

5. **Analyze & package lightweight results to send back:**
   ```bash
   python 03_analyze_ucst.py
   bash 04_package_results.sh
   ```
   Send back **`PNAGA40_UCST_RESULTS.tar.gz`** (typically ~50–100 MB; excludes huge full-water trajectories).
