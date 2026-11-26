import glob
import os
import shutil
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = os.path.abspath(os.getcwd())
WORK = os.path.join(ROOT, "_build_tmp")
OUT_DIR = os.path.join(ROOT, "00_master_inputs")
if os.path.exists(WORK):
    shutil.rmtree(WORK)
os.makedirs(WORK, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)
os.chdir(WORK)

# Auto-detect gmx or gmx_mpi
gmx_bin = os.environ.get("GMX_BIN")
if not gmx_bin:
    gmx_bin = shutil.which("gmx") or shutil.which("gmx_mpi")
if not gmx_bin:
    raise RuntimeError("Neither 'gmx' nor 'gmx_mpi' found in PATH. Load your GROMACS module first.")
print(f"Using GROMACS binary: {gmx_bin}")

# 1. Build methyl-capped NAGA repeat unit (exact double bonds preserved via SDF)
smiles = "CCC(C)C(=O)NCC(=O)N"
mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(mol)
Chem.MolToMolFile(mol, "capped_naga.sdf")

pat = Chem.MolFromSmarts("[CH3][CH2][CH]([CH3])[C](=[O])[NH][CH2][C](=[O])[NH2]")
match = mol.GetSubstructMatch(pat)
assert len(match) == 11, "Atom mapping failed!"

c1_idx, c2_idx, c3_idx, c4_idx = match[0], match[1], match[2], match[3]
left_cap_ids = {c1_idx + 1} | {
    n.GetIdx() + 1
    for n in mol.GetAtomWithIdx(c1_idx).GetNeighbors()
    if n.GetAtomicNum() == 1
}
right_cap_ids = {c4_idx + 1} | {
    n.GetIdx() + 1
    for n in mol.GetAtomWithIdx(c4_idx).GetNeighbors()
    if n.GetAtomicNum() == 1
}
c2_id, c3_id = c2_idx + 1, c3_idx + 1


def run(cmd):
    print(f">>> {cmd}")
    subprocess.run(cmd, shell=True, check=True)


# 2. Generate exact GAFF2 atom types + AM1-BCC partial charges directly to MOL2
run("antechamber -i capped_naga.sdf -fi mdl -o capped_naga.mol2 -fo mol2 -c bcc -at gaff2 -nc 0 -rn MNG")
run("parmchk2 -i capped_naga.mol2 -f mol2 -o naga.frcmod -s gaff2")


# 3. Slice into Head (HNG), Middle (MNG), and Tail (TNG) MOL2 blocks
def parse_mol2(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    atoms, bonds, mode = [], [], None
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            mode = "ATOM"
            continue
        if line.startswith("@<TRIPOS>BOND"):
            mode = "BOND"
            continue
        if line.startswith("@<TRIPOS>"):
            mode = None
            continue
        if not line.strip():
            continue
        p = line.split()
        if mode == "ATOM":
            atoms.append(
                {
                    "id": int(p[0]),
                    "name": p[1],
                    "x": float(p[2]),
                    "y": float(p[3]),
                    "z": float(p[4]),
                    "type": p[5],
                    "q": float(p[8]),
                }
            )
        elif mode == "BOND":
            bonds.append((int(p[1]), int(p[2]), p[3]))
    return atoms, bonds


atoms, bonds = parse_mol2("capped_naga.mol2")
print("GAFF2 atom types assigned:", [a["type"] for a in atoms[:11]])


def write_residue_mol2(filename, resname, omit_ids):
    kept = [dict(a) for a in atoms if a["id"] not in omit_ids]
    id_map = {a["id"]: i + 1 for i, a in enumerate(kept)}
    for a in kept:
        if a["id"] == c2_id:
            a["name"] = "C2"
        if a["id"] == c3_id:
            a["name"] = "C3"
    shift = -sum(a["q"] for a in kept) / len(kept)
    for a in kept:
        a["q"] = round(a["q"] + shift, 6)
    kept[-1]["q"] = round(kept[-1]["q"] - sum(a["q"] for a in kept), 6)

    kept_bonds = [
        (id_map[a1], id_map[a2], bt)
        for (a1, a2, bt) in bonds
        if a1 in id_map and a2 in id_map
    ]
    with open(filename, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{resname}\n")
        f.write(f" {len(kept)} {len(kept_bonds)} 1 0 0\nSMALL\nbcc\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for i, a in enumerate(kept, 1):
            f.write(
                f"{i:7d} {a['name']:<8s} {a['x']:10.4f} {a['y']:10.4f} {a['z']:10.4f} "
                f"{a['type']:<6s} 1 {resname:<4s} {a['q']:10.6f}\n"
            )
        f.write("@<TRIPOS>BOND\n")
        for i, (b1, b2, bt) in enumerate(kept_bonds, 1):
            f.write(f"{i:6d} {b1:5d} {b2:5d} {bt}\n")
        f.write(f"@<TRIPOS>SUBSTRUCTURE\n     1 {resname} 1 TEMP 0 **** **** 0 ROOT\n")


write_residue_mol2("HNG.mol2", "HNG", right_cap_ids)
write_residue_mol2("MNG.mol2", "MNG", left_cap_ids | right_cap_ids)
write_residue_mol2("TNG.mol2", "TNG", left_cap_ids)

# 4. Assemble 40-mer PNAGA in tleap (1 Head + 38 Middle + 1 Tail)
seq = " ".join(["HNG"] + ["MNG"] * 38 + ["TNG"])
with open("leap.in", "w", encoding="utf-8", newline="\n") as f:
    f.write(
        f"""source leaprc.gaff2
loadamberparams naga.frcmod
HNG = loadmol2 HNG.mol2
MNG = loadmol2 MNG.mol2
TNG = loadmol2 TNG.mol2
set HNG tail HNG.1.C3
set MNG head MNG.1.C2
set MNG tail MNG.1.C3
set TNG head TNG.1.C2
poly = sequence {{ {seq} }}
check poly
savepdb poly pnaga40.pdb
saveamberparm poly pnaga40.prmtop pnaga40.inpcrd
quit
"""
    )
run("tleap -f leap.in")

# 5. Convert AMBER topology to GROMACS
run("acpype -p pnaga40.prmtop -x pnaga40.inpcrd -b pnaga40")

gro_file = glob.glob("**/*GMX.gro", recursive=True)[0]
top_file = glob.glob("**/*GMX.top", recursive=True)[0]
shutil.copy(gro_file, os.path.join(OUT_DIR, "pnaga40_GMX.gro"))
shutil.copy("pnaga40.pdb", os.path.join(OUT_DIR, "pnaga40.pdb"))

# 6. Inject AMBER/GAFF2 TIP3P water parameters into topol_master.top and solvate
top_text = open(top_file, encoding="utf-8").read()
tip3p_atomtypes = (
    " OW          OW          0.00000  0.00000   A     3.15061e-01   6.36386e-01 ; TIP3P O\n"
    " HW          HW          0.00000  0.00000   A     0.00000e+00   0.00000e+00 ; TIP3P H\n"
)
tip3p_mol = """
[ moleculetype ]
; molname       nrexcl
SOL             2

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   OW          1       SOL       OW       1     -0.834    16.00000
  2   HW          1       SOL       HW1      1      0.417     1.00800
  3   HW          1       SOL       HW2      1      0.417     1.00800

[ settles ]
; i     funct   doh        dhh
1       1       0.09572    0.15139

[ exclusions ]
1 2 3
2 1 3
3 1 2

"""
if "TIP3P O" not in top_text:
    top_text = top_text.replace("[ moleculetype ]", tip3p_atomtypes + "\n[ moleculetype ]", 1)
    top_text = top_text.replace("[ system ]", tip3p_mol + "[ system ]", 1)

with open(os.path.join(OUT_DIR, "topol_master.top"), "w", encoding="utf-8", newline="\n") as f:
    f.write(top_text)

os.chdir(OUT_DIR)
run(f'"{gmx_bin}" editconf -f pnaga40_GMX.gro -o box.gro -c -d 1.5 -bt dodecahedron')
run(f'"{gmx_bin}" solvate -cp box.gro -cs spc216.gro -o solv_master.gro -p topol_master.top')
shutil.rmtree(WORK, ignore_errors=True)
print("SUCCESS! Master solvated 40-mer ready in 00_master_inputs/")
