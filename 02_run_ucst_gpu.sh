#!/usr/bin/env bash
set -euo pipefail

# Target production length in nanoseconds (change to 200 if desired)
TARGET_NS=100
NSTEPS=$((TARGET_NS * 500000))  # dt = 0.002 ps

# Pass one temperature as argument (e.g. "bash 02_run_ucst_gpu.sh 30") or leave blank for all 6
if [ "$#" -ge 1 ]; then
    TEMPS=("$@")
else
    TEMPS=(30 40 50 60 70 80)
fi

ROOT_DIR="$(pwd)"
MASTER_DIR="${ROOT_DIR}/00_master_inputs"
GMX_BIN="${GMX_BIN:-gmx}"

if [ ! -f "${MASTER_DIR}/solv_master.gro" ]; then
    echo "ERROR: ${MASTER_DIR}/solv_master.gro not found. Run 'python 01_build_pnaga40.py' first."
    exit 1
fi

for TEMP_C in "${TEMPS[@]}"; do
    TEMP_K=$(awk "BEGIN {printf \"%.2f\", ${TEMP_C} + 273.15}")
    RUN_DIR="${ROOT_DIR}/runs/T${TEMP_C}C"
    mkdir -p "${RUN_DIR}"
    cd "${RUN_DIR}"
    echo "============================================================"
    echo "Running 40-mer PNAGA at ${TEMP_C} C (${TEMP_K} K) for ${TARGET_NS} ns"
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
nstxout-compressed = 10000 ; save every 20 ps
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

    if [ ! -f em.gro ]; then
        ${GMX_BIN} grompp -f em.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 1
        ${GMX_BIN} mdrun -deffnm em
    fi

    if [ ! -f nvt.gro ]; then
        ${GMX_BIN} grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 1
        ${GMX_BIN} mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu -update gpu
    fi

    if [ ! -f npt.gro ]; then
        ${GMX_BIN} grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
        ${GMX_BIN} mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu -update gpu
    fi

    if [ ! -f prod.tpr ]; then
        ${GMX_BIN} grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -o prod.tpr -maxwarn 1
    fi

    # Production MD with auto-resume checkpointing every 15 min
    ${GMX_BIN} mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu -update gpu -cpi prod.cpt -cpt 15

    # Post-processing: whole polymer trajectory, Rg, SASA, and H-bonds
    echo "non-Water" | ${GMX_BIN} trjconv -s prod.tpr -f prod.xtc -o prod_polymer.xtc -pbc whole
    echo "non-Water" | ${GMX_BIN} gyrate -s prod.tpr -f prod_polymer.xtc -o gyrate.xvg
    echo "non-Water" | ${GMX_BIN} sasa -s prod.tpr -f prod_polymer.xtc -o sasa.xvg
    if ${GMX_BIN} help hbond-legacy >/dev/null 2>&1; then
        printf "non-Water\nnon-Water\n" | ${GMX_BIN} hbond-legacy -s prod.tpr -f prod_polymer.xtc -num hbond_intra.xvg
    else
        printf "non-Water\nnon-Water\n" | ${GMX_BIN} hbond -s prod.tpr -f prod_polymer.xtc -num hbond_intra.xvg || true
    fi

    cd "${ROOT_DIR}"
done

echo "ALL REQUESTED TEMPERATURES COMPLETED."
