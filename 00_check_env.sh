#!/usr/bin/env bash
set -euo pipefail

echo "=== 1. GPU Check ==="
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,memory.total --format=csv
else
    echo "WARNING: nvidia-smi not found on this login node (normal if GPUs are only on compute nodes)."
fi

echo ""
echo "=== 2. GROMACS Check ==="
if command -v gmx >/dev/null 2>&1; then
    GMX_BIN="gmx"
elif command -v gmx_mpi >/dev/null 2>&1; then
    GMX_BIN="gmx_mpi"
else
    echo "ERROR: Neither 'gmx' nor 'gmx_mpi' found in PATH. Load your GROMACS module first."
    exit 1
fi
echo "Found GROMACS binary: ${GMX_BIN}"
${GMX_BIN} --version | grep -E "GROMACS version|GPU support|CUDA" || true

echo ""
echo "=== 3. Builder & Analysis Tools Check ==="
for tool in python antechamber parmchk2 tleap acpype; do
    if command -v "${tool}" >/dev/null 2>&1; then
        echo "  [OK] ${tool}: $(command -v ${tool})"
    else
        echo "  [MISSING] ${tool} (run: conda activate pnaga_md)"
    fi
done

python -c "import rdkit, numpy, pandas, matplotlib; print('  [OK] Python packages: rdkit, numpy, pandas, matplotlib')"
echo ""
echo "Pre-flight check complete."
