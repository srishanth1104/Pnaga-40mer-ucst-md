#!/usr/bin/env bash
set -euo pipefail

OUT="PNAGA40_UCST_RESULTS.tar.gz"
ITEMS=()

[ -d "analysis" ] && ITEMS+=("analysis")
[ -f "00_master_inputs/pnaga40.pdb" ] && ITEMS+=("00_master_inputs/pnaga40.pdb")

for rdir in runs/T*C; do
    [ -d "${rdir}" ] || continue
    for f in gyrate.xvg sasa.xvg hbond_intra.xvg hbond_water.xvg npt_temp.xvg npt_press.xvg npt_dens.xvg em_pot.xvg polymer.gro polymer.tpr prod_polymer.xtc prod.log; do
        [ -f "${rdir}/${f}" ] && ITEMS+=("${rdir}/${f}")
    done
done

if [ "${#ITEMS[@]}" -eq 0 ]; then
    echo "ERROR: No result files found to package."
    exit 1
fi

tar -czvf "${OUT}" "${ITEMS[@]}"
echo ""
echo "READY TO SEND: ${OUT} ($(du -h "${OUT}" | cut -f1))"
