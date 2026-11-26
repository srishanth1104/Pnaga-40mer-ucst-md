#!/usr/bin/env bash
set -euo pipefail

# === Configurable Parameters ===
TARGET_NS="${TARGET_NS:-100}"     # Production length per temperature in ns (default 100 ns)
MAX_HOURS="${MAX_HOURS:-47.5}"    # Clean checkpoint exit before cluster walltime (default 47.5 h)
NSTEPS=$((TARGET_NS * 500000))    # dt = 0.002 ps

# Auto-detect gmx or gmx_mpi
if [ -n "${GMX_BIN:-}" ]; then
    GMX="${GMX_BIN}"
elif command -v gmx >/dev/null 2>&1; then
    GMX="gmx"
elif command -v gmx_mpi >/dev/null 2>&1; then
    GMX="gmx_mpi"
else
    echo "ERROR: Neither 'gmx' nor 'gmx_mpi' found in PATH. Load your GROMACS GPU module first."
    exit 1
fi

# Temperatures: pass one or more as arguments (e.g. "bash 02_run_ucst_gpu.sh 30") or leave blank for all 6
if [ "$#" -ge 1 ]; then
    TEMPS=("$@")
else
    TEMPS=(30 40 50 60 70 80)
fi

ROOT_DIR="$(pwd)"
MASTER_DIR="${ROOT_DIR}/00_master_inputs"

if [ ! -f "${MASTER_DIR}/solv_master.gro" ]; then
    echo "ERROR: ${MASTER_DIR}/solv_master.gro not found. Run 'python 01_build_pnaga40.py' first."
    exit 1
fi

# Helper: try full GPU offload first, fall back cleanly to -nb gpu -pme gpu if older GROMACS rejects -update gpu
run_gpu_mdrun() {
    local deffnm="$1"
    shift
    if ! ${GMX} mdrun -deffnm "${deffnm}" -nb gpu -pme gpu -bonded gpu -update gpu "$@"; then
        echo "Falling back to -nb gpu -pme gpu for ${deffnm}..."
        ${GMX} mdrun -deffnm "${deffnm}" -nb gpu -pme gpu "$@"
    fi
}

for TEMP_C in "${TEMPS[@]}"; do
    TEMP_K=$(awk "BEGIN {printf \"%.2f\", ${TEMP_C} + 273.15}")
    RUN_DIR="${ROOT_DIR}/runs/T${TEMP_C}C"
    mkdir -p "${RUN_DIR}"
    cd "${RUN_DIR}"
    echo "============================================================"
    echo "Running 40-mer PNAGA at ${TEMP_C} C (${TEMP_K} K) for ${TARGET_NS} ns using ${GMX}"
    echo "============================================================"

    cp -f "${MASTER_DIR}/solv_master.gro" solv.gro
    cp -f "${MASTER_DIR}/topol_master.top" topol.top

    cat > em.mdp << 'EOF'
integrator    = steep
emtol         = 1000.0
emstep        = 0.01
nsteps        = 50000
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
pbc           = xyz
EOF

    cat > nvt.mdp << EOF
integrator    = md
dt            = 0.002
nsteps        = 100000     ; 200 ps
nstxout-compressed = 5000
nstenergy     = 1000
nstlog        = 1000
continuation  = no
constraint_algorithm = lincs
constraints   = h-bonds
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
DispCorr      = EnerPres
tcoupl        = V-rescale
tc-grps       = System
tau_t         = 0.5
ref_t         = ${TEMP_K}
pcoupl        = no
pbc           = xyz
gen_vel       = yes
gen_temp      = ${TEMP_K}
gen_seed      = -1
EOF

    cat > npt.mdp << EOF
integrator    = md
dt            = 0.002
nsteps        = 500000     ; 1 ns
nstxout-compressed = 5000
nstenergy     = 2500
nstlog        = 2500
continuation  = yes
constraint_algorithm = lincs
constraints   = h-bonds
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
DispCorr      = EnerPres
tcoupl        = V-rescale
tc-grps       = System
tau_t         = 0.5
ref_t         = ${TEMP_K}
pcoupl        = C-rescale
pcoupltype    = isotropic
tau_p         = 2.0
ref_p         = 1.0
compressibility = 4.5e-5
pbc           = xyz
gen_vel       = no
EOF

    cat > prod.mdp << EOF
integrator    = md
dt            = 0.002
nsteps        = ${NSTEPS}
nstxout-compressed = 10000 ; save frame every 20 ps
nstenergy     = 5000
nstlog        = 5000
continuation  = yes
constraint_algorithm = lincs
constraints   = h-bonds
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
DispCorr      = EnerPres
tcoupl        = V-rescale
tc-grps       = System
tau_t         = 0.5
ref_t         = ${TEMP_K}
pcoupl        = C-rescale
pcoupltype    = isotropic
tau_p         = 2.0
ref_p         = 1.0
compressibility = 4.5e-5
pbc           = xyz
gen_vel       = no
EOF

    # 1. Energy Minimization
    if [ ! -f em.gro ]; then
        ${GMX} grompp -f em.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 1
        ${GMX} mdrun -deffnm em
        printf "Potential\n0\n" | ${GMX} energy -f em.edr -o em_pot.xvg || true
    fi

    # 2. NVT Equilibration (200 ps)
    if [ ! -f nvt.gro ]; then
        ${GMX} grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 1
        run_gpu_mdrun nvt
    fi

    # 3. NPT Equilibration (1 ns) + QC extraction
    if [ ! -f npt.gro ]; then
        ${GMX} grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
        run_gpu_mdrun npt
        printf "Temperature\n0\n" | ${GMX} energy -f npt.edr -o npt_temp.xvg || true
        printf "Pressure\n0\n" | ${GMX} energy -f npt.edr -o npt_press.xvg || true
        printf "Density\n0\n" | ${GMX} energy -f npt.edr -o npt_dens.xvg || true
    fi

    # 4. Production MD (auto-resumes seamlessly via prod.cpt)
    if [ ! -f prod.tpr ]; then
        ${GMX} grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -o prod.tpr -maxwarn 1
    fi

    run_gpu_mdrun prod -cpi prod.cpt -cpt 15 -maxh "${MAX_HOURS}"

    # 5. Post-processing (only runs once prod.gro exists, meaning production completed)
    if [ -f prod.gro ]; then
        echo "Extracting polymer-only trajectory and UCST observables..."
        # Create deterministic index with explicit Polymer (! r SOL) and Water (r SOL) groups
        printf "! r SOL\nname 10 Polymer\nq\n" | ${GMX} make_ndx -f prod.tpr -o index.ndx || true

        # Polymer-only matching tpr, gro, and whole-molecule xtc (zero atom-count mismatch!)
        echo "non-Water" | ${GMX} convert-tpr -s prod.tpr -o polymer.tpr
        echo "non-Water" | ${GMX} trjconv -s prod.tpr -f prod.gro -o polymer.gro -pbc whole
        echo "non-Water" | ${GMX} trjconv -s prod.tpr -f prod.xtc -o prod_polymer.xtc -pbc whole

        # Rg and SASA on matching polymer.tpr + prod_polymer.xtc
        echo "System" | ${GMX} gyrate -s polymer.tpr -f prod_polymer.xtc -o gyrate.xvg
        echo "System" | ${GMX} sasa -s polymer.tpr -f prod_polymer.xtc -o sasa.xvg

        # Intra-chain (polymer-polymer) and polymer-water hydrogen bonds
        if ${GMX} help hbond-legacy >/dev/null 2>&1; then
            printf "System\nSystem\n" | ${GMX} hbond-legacy -s polymer.tpr -f prod_polymer.xtc -num hbond_intra.xvg
            printf "non-Water\nWater\n" | ${GMX} hbond-legacy -s prod.tpr -f prod.xtc -num hbond_water.xvg || true
        else
            printf "System\nSystem\n" | ${GMX} hbond -s polymer.tpr -f prod_polymer.xtc -num hbond_intra.xvg || true
            printf "non-Water\nWater\n" | ${GMX} hbond -s prod.tpr -f prod.xtc -num hbond_water.xvg || true
        fi
    else
        echo "NOTE: prod.gro not yet written (stopped cleanly at MAX_HOURS). Re-submit to auto-resume from prod.cpt."
    fi

    cd "${ROOT_DIR}"
done

echo "DONE."
