# 40-mer PNAGA UCST Molecular Dynamics (GROMACS GPU)

End-to-end GROMACS GPU workflow to study the Upper Critical Solution Temperature (UCST) coil-to-globule transition of a single 40-mer poly(*N*-acryloyl glycinamide) (PNAGA) chain in explicit TIP3P water across **30, 40, 50, 60, 70, and 80 °C** (100 ns production per temperature, GAFF2 + AM1-BCC charges).

## Quick Start on a GPU Cluster

1. **Create the builder/analysis environment (one-time):**
   ```bash
   conda env create -f environment.yml
   conda activate pnaga_md
   ```
   *(Load your cluster's GPU-enabled GROMACS module, e.g., `module load gromacs`, before running Steps 2 and 3.)*

2. **Build and solvate the master 40-mer (~1 minute):**
   ```bash
   python 01_build_pnaga40.py
   ```
   Generates `00_master_inputs/solv_master.gro` and `00_master_inputs/topol_master.top` (40-mer PNAGA in a rhombic dodecahedron TIP3P water box with 1.5 nm padding).

3. **Run all 6 temperatures on GPU:**
   - **Option A (SLURM GPU array — runs all 6 temperatures in parallel):** Edit `#SBATCH --partition` and `module load` in `submit_slurm.sh`, then submit:
     ```bash
     sbatch submit_slurm.sh
     ```
   - **Option B (Single GPU sequential run for all 6 temperatures):**
     ```bash
     bash 02_run_ucst_gpu.sh
     ```
   - **Option C (Run one specific temperature, e.g. 30 °C):**
     ```bash
     bash 02_run_ucst_gpu.sh 30
     ```

4. **Plot the UCST transition curve once finished:**
   ```bash
   python 03_analyze_ucst.py
   ```
   Outputs `analysis/ucst_summary.csv` and `analysis/pnaga40_ucst_curve.png`.
