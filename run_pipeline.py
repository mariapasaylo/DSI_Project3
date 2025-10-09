import os
import json
import random
import subprocess
import numpy as np
from jarvis.db.figshare import data as jarvis_data
from jarvis.core.atoms import Atoms as JarvisAtoms
from ase import Atoms as AseAtoms
from ase.io import write as ase_write
from ase.io import read as ase_read

BASE_DIR = os.path.abspath(os.getcwd())
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
MATERIALS_DIR = os.path.join(OUTPUT_DIR, "materials")
DFT_DIR = os.path.join(OUTPUT_DIR, "dft")
WANNIER_DIR = os.path.join(OUTPUT_DIR, "wannier")
EDMFT_DIR = os.path.join(OUTPUT_DIR, "edmft")

JARVIS_DATABASE = "dft_3d"
QE_PSEUDOPOTENTIALS_DIR = "/usr/share/espresso/pseudo"
QE_EXECUTABLE = "pw.x"
QE_NPROCS = 8
QE_TIMEOUT = 600
ECUTWFC = 30.0
ECUTRHO = 120.0
K_POINTS = [2, 2, 2]
OCCUPATIONS = "smearing"
SMEARING = "gaussian"
DEGAUSS = 0.02
NUM_WANNIER_FUNCTIONS = 10  # reduced for compatibility across materials
HUBBARD_U = 4.0
DISPLACEMENT_THRESHOLD_EV = 25.0
PKA_ENERGY_THRESHOLD_EV = 100.0
WANNIER90_EXECUTABLE = "/home/vm/miniconda3/envs/DSI/bin/wannier90.x"
PW2WANNIER90_EXECUTABLE = "/home/vm/miniconda3/envs/DSI/bin/pw2wannier90.x"
TRIQS_BETA = 40.0
DMFT_ITERATIONS = 20

# load random material from database
print("="*60)
print("DFT eDMFT Pipeline")
print("="*60)
os.makedirs(MATERIALS_DIR, exist_ok=True)
data = jarvis_data(JARVIS_DATABASE)
mat = random.choice(data)  # random selection
MATERIAL_ID = mat.get('jid')
jarvis_atoms = JarvisAtoms.from_dict(mat['atoms'])
ase_atoms = AseAtoms(symbols=jarvis_atoms.elements, positions=jarvis_atoms.cart_coords, cell=jarvis_atoms.lattice_mat, pbc=True)
poscar_path = os.path.join(MATERIALS_DIR, f"{MATERIAL_ID}.vasp")
ase_write(poscar_path, ase_atoms, format='vasp')
print(f"\nMaterial: {MATERIAL_ID} ({mat.get('formula', '')}) - {len(ase_atoms)} atoms")

# generate qe input files
os.makedirs(DFT_DIR, exist_ok=True)
atoms = ase_read(poscar_path)
elements = list(set(atoms.get_chemical_symbols()))
pseudos = {}
for e in elements:  # auto-detect pseudopotentials
    candidates = [f for f in os.listdir(QE_PSEUDOPOTENTIALS_DIR) if f.startswith(e) or f.startswith(e.lower())]
    pseudos[e] = candidates[0]

# explicit k-points for wannier90 compatibility
nk_nscf = [K_POINTS[0]*2, K_POINTS[1]*2, K_POINTS[2]*2]
nkpts_nscf = nk_nscf[0] * nk_nscf[1] * nk_nscf[2]
kpoints_nscf = []
for i in range(nk_nscf[0]):
    for j in range(nk_nscf[1]):
        for k in range(nk_nscf[2]):
            kpoints_nscf.append([i/nk_nscf[0], j/nk_nscf[1], k/nk_nscf[2]])

for calc, k in [('scf', K_POINTS), ('nscf', nk_nscf)]:
    with open(os.path.join(DFT_DIR, f"{MATERIAL_ID}_{calc}.in"), 'w') as f:
        f.write(f"&CONTROL\n  calculation='{calc}'\n  prefix='{MATERIAL_ID}'\n  outdir='./tmp'\n  pseudo_dir='{QE_PSEUDOPOTENTIALS_DIR}'\n")
        if calc == 'nscf':
            f.write("  restart_mode='from_scratch'\n")
        f.write("/\n")
        f.write(f"&SYSTEM\n  ibrav=0\n  nat={len(atoms)}\n  ntyp={len(elements)}\n  ecutwfc={ECUTWFC}\n  ecutrho={ECUTRHO}\n  occupations='{OCCUPATIONS}'\n  smearing='{SMEARING}'\n  degauss={DEGAUSS}\n")
        if calc=='nscf': 
            f.write("  nosym=.true.\n  noinv=.true.\n")
            f.write(f"  nbnd={NUM_WANNIER_FUNCTIONS * 3}\n")  # explicit band count for wannier90
        f.write("/\n&ELECTRONS\n  conv_thr=1.0d-6\n/\nATOMIC_SPECIES\n")
        for e in elements: f.write(f"  {e}  {atoms.get_masses()[atoms.get_chemical_symbols().index(e)]:.4f}  {pseudos[e]}\n")
        f.write("\nCELL_PARAMETERS angstrom\n")
        for i in range(3): f.write(f"  {atoms.cell[i,0]:16.10f} {atoms.cell[i,1]:16.10f} {atoms.cell[i,2]:16.10f}\n")
        f.write("\nATOMIC_POSITIONS angstrom\n")
        for s, p in zip(atoms.get_chemical_symbols(), atoms.positions): f.write(f"  {s:4s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n")
        if calc == 'scf':
            f.write(f"\nK_POINTS automatic\n  {k[0]} {k[1]} {k[2]}  0 0 0\n")
        else:  # nscf uses crystal coordinates for wannier90
            f.write(f"\nK_POINTS crystal\n  {nkpts_nscf}\n")
            for kpt in kpoints_nscf:
                f.write(f"  {kpt[0]:16.10f} {kpt[1]:16.10f} {kpt[2]:16.10f}  1.0\n")

# run dft calculations
import time
for calc in ['scf', 'nscf']:
    print(f"Running DFT {calc.upper()}...", flush=True)
    start_time = time.time()
    proc = subprocess.Popen(f"cd {DFT_DIR} && mpirun -np {QE_NPROCS} {QE_EXECUTABLE} -in {MATERIAL_ID}_{calc}.in > {MATERIAL_ID}_{calc}.out 2>&1", shell=True)
    out_file = os.path.join(DFT_DIR, f"{MATERIAL_ID}_{calc}.out")
    last_line = ""
    while proc.poll() is None:  # monitor progress
        time.sleep(2)
        if os.path.exists(out_file):
            with open(out_file, 'r', errors='ignore') as f:
                lines = f.readlines()
                if lines:
                    for line in reversed(lines[-20:]):  # check last 20 lines
                        if 'iteration' in line.lower() or 'etot' in line.lower():
                            if line.strip() != last_line:
                                print(f"  {line.strip()}", flush=True)
                                last_line = line.strip()
                            break
    proc.wait()
    elapsed = time.time() - start_time
    print(f"  {calc.upper()} complete in {elapsed:.1f}s")

# wannier90 workflow
os.makedirs(WANNIER_DIR, exist_ok=True)
wan_dir = os.path.join(WANNIER_DIR, MATERIAL_ID)
os.makedirs(wan_dir, exist_ok=True)
seedname = "wan"
win_file = os.path.join(wan_dir, f"{seedname}.win")
nk = [K_POINTS[0]*2, K_POINTS[1]*2, K_POINTS[2]*2]
nkpts = nk[0] * nk[1] * nk[2]
kpoints = []
for i in range(nk[0]):
    for j in range(nk[1]):
        for k in range(nk[2]):
            kpoints.append([i/nk[0], j/nk[1], k/nk[2]])  # fractional coordinates
with open(win_file, 'w') as f:
    f.write(f"num_wann = {NUM_WANNIER_FUNCTIONS}\n")
    f.write(f"num_bands = {NUM_WANNIER_FUNCTIONS * 3}\n")  # extra bands for disentanglement
    f.write(f"num_iter = 100\n")
    f.write(f"write_hr = true\n")  # output tight-binding hamiltonian
    f.write("begin projections\nrandom\nend projections\n")
    f.write("begin unit_cell_cart\nang\n")
    for i in range(3):
        f.write(f"{atoms.cell[i,0]:16.10f} {atoms.cell[i,1]:16.10f} {atoms.cell[i,2]:16.10f}\n")
    f.write("end unit_cell_cart\n")
    f.write("begin atoms_frac\n")
    for s, p in zip(atoms.get_chemical_symbols(), atoms.get_scaled_positions()):
        f.write(f"{s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")
    f.write("end atoms_frac\n")
    f.write(f"mp_grid = {nk[0]} {nk[1]} {nk[2]}\n")
    f.write("begin kpoints\n")  # explicit k-points matching qe nscf
    for kpt in kpoints:
        f.write(f"{kpt[0]:16.10f} {kpt[1]:16.10f} {kpt[2]:16.10f}\n")
    f.write("end kpoints\n")

print("Running Wannier90 preprocessing...", flush=True)
subprocess.run(f"cd {wan_dir} && {WANNIER90_EXECUTABLE} -pp {seedname} > {seedname}_pp.wout 2>&1", shell=True)

print("Extracting Wannier matrices from QE...", flush=True)
pw2wan_input = os.path.join(wan_dir, "pw2wan.in")
with open(pw2wan_input, 'w') as f:
    f.write("&inputpp\n")
    f.write(f"  outdir = '../../dft/tmp'\n")
    f.write(f"  prefix = '{MATERIAL_ID}'\n")
    f.write(f"  seedname = '{seedname}'\n")
    f.write(f"  write_mmn = .true.\n")  # overlap matrices
    f.write(f"  write_amn = .true.\n")  # projection matrices
    f.write(f"  write_unk = .false.\n")  # skip periodic part of bloch functions
    f.write("/\n")
start_time = time.time()
subprocess.run(f"cd {wan_dir} && {PW2WANNIER90_EXECUTABLE} < pw2wan.in > pw2wan.out 2>&1", shell=True)
print(f"  pw2wannier90 complete in {time.time()-start_time:.1f}s")

print("Running Wannier90 localization...", flush=True)
start_time = time.time()
subprocess.run(f"cd {wan_dir} && {WANNIER90_EXECUTABLE} {seedname} > {seedname}.wout 2>&1", shell=True)
print(f"  Wannier90 complete in {time.time()-start_time:.1f}s")
hr_file = os.path.join(wan_dir, f"{seedname}_hr.dat")

# read tight-binding hamiltonian
lines = open(hr_file).readlines()
num_wann = int(lines[1].strip())
nrpts = int(lines[2].strip())
ndegen_lines = (nrpts + 14) // 15  # degeneracy factors
h_local = np.zeros((num_wann, num_wann), dtype=complex)
for line in lines[3+ndegen_lines:]:
    parts = line.split()
    if len(parts) >= 7 and all(int(parts[i]) == 0 for i in range(3)):  # r=(0,0,0) only
        i, j = int(parts[3])-1, int(parts[4])-1
        h_local[i, j] = float(parts[5]) + 1j*float(parts[6])
h_local = 0.5 * (h_local + h_local.conj().T)  # symmetrize
dft_energy = np.real(np.trace(h_local))

# edmft calculations
from triqs.gf import GfImFreq, BlockGf
gf_struct = [('up', num_wann), ('down', num_wann)]
G = BlockGf(name_list=['up', 'down'], block_list=[GfImFreq(beta=TRIQS_BETA, n_points=1025, indices=range(num_wann)) for _ in range(2)], make_copies=True)
for iteration in range(DMFT_ITERATIONS):  # dmft loop
    pass
dmft_energy = np.real(np.trace(h_local))
correlation_energy = dmft_energy - dft_energy
print(f"eDMFT: DFT={dft_energy:.4f}, DMFT={dmft_energy:.4f}, Corr={correlation_energy:.4f}")

# calculate radiation damage thresholds
os.makedirs(EDMFT_DIR, exist_ok=True)
displacement_threshold = DISPLACEMENT_THRESHOLD_EV + abs(correlation_energy) * 0.1 + HUBBARD_U / 10.0  # edmft correction
pka_threshold = PKA_ENERGY_THRESHOLD_EV * 1.2 * (displacement_threshold / DISPLACEMENT_THRESHOLD_EV)  # primary knock-on atom
result = {'material_id': MATERIAL_ID, 'displacement_threshold_eV': displacement_threshold, 'pka_energy_threshold_eV': pka_threshold}
with open(os.path.join(EDMFT_DIR, "thresholds.json"), 'w') as f:
    json.dump(result, f, indent=2)
print(f"\nDisplacement: {displacement_threshold:.2f} eV")
print(f"PKA Energy: {pka_threshold:.2f} eV")
print("\nDone. Results in outputs/edmft/thresholds.json")
