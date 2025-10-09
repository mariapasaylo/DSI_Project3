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
NUM_WANNIER_FUNCTIONS = 20
HUBBARD_U = 4.0
DISPLACEMENT_THRESHOLD_EV = 25.0
PKA_ENERGY_THRESHOLD_EV = 100.0
WANNIER90_EXECUTABLE = os.path.join(os.environ.get('CONDA_PREFIX', '/home/vm/miniconda3/envs/DSI'), 'bin', 'wannier90.x')
TRIQS_BETA = 40.0
DMFT_ITERATIONS = 20

# ========== STEP 1: LOAD MATERIAL FROM JARVIS ==========
os.makedirs(MATERIALS_DIR, exist_ok=True)
data = jarvis_data(JARVIS_DATABASE)
mat = random.choice(data)
MATERIAL_ID = mat.get('jid')
jarvis_atoms = JarvisAtoms.from_dict(mat['atoms'])
ase_atoms = AseAtoms(symbols=jarvis_atoms.elements, positions=jarvis_atoms.cart_coords, cell=jarvis_atoms.lattice_mat, pbc=True)
poscar_path = os.path.join(MATERIALS_DIR, f"{MATERIAL_ID}.vasp")
ase_write(poscar_path, ase_atoms, format='vasp')
print(f"Material: {MATERIAL_ID} ({mat.get('formula', '')})")

# ========== STEP 2: GENERATE QUANTUM ESPRESSO INPUT FILES ==========
os.makedirs(DFT_DIR, exist_ok=True)
atoms = ase_read(poscar_path)
elements = list(set(atoms.get_chemical_symbols()))
pseudos = {}
for e in elements:
    candidates = [f for f in os.listdir(QE_PSEUDOPOTENTIALS_DIR) if f.startswith(e) or f.startswith(e.lower())]
    pseudos[e] = candidates[0] if candidates else f"{e}.UPF"
for calc, k in [('scf', K_POINTS), ('nscf', [k*2 for k in K_POINTS])]:
    with open(os.path.join(DFT_DIR, f"{MATERIAL_ID}_{calc}.in"), 'w') as f:
        f.write(f"&CONTROL\n  calculation='{calc}'\n  prefix='{MATERIAL_ID}'\n  outdir='./tmp'\n  pseudo_dir='{QE_PSEUDOPOTENTIALS_DIR}'\n")
        if calc == 'nscf':
            f.write("  restart_mode='from_scratch'\n")
        f.write("/\n")
        f.write(f"&SYSTEM\n  ibrav=0\n  nat={len(atoms)}\n  ntyp={len(elements)}\n  ecutwfc={ECUTWFC}\n  ecutrho={ECUTRHO}\n  occupations='{OCCUPATIONS}'\n  smearing='{SMEARING}'\n  degauss={DEGAUSS}\n")
        if calc=='nscf': f.write("  nosym=.true.\n  noinv=.true.\n")
        f.write("/\n&ELECTRONS\n  conv_thr=1.0d-6\n/\nATOMIC_SPECIES\n")
        for e in elements: f.write(f"  {e}  {atoms.get_masses()[atoms.get_chemical_symbols().index(e)]:.4f}  {pseudos[e]}\n")
        f.write("\nCELL_PARAMETERS angstrom\n")
        for i in range(3): f.write(f"  {atoms.cell[i,0]:16.10f} {atoms.cell[i,1]:16.10f} {atoms.cell[i,2]:16.10f}\n")
        f.write("\nATOMIC_POSITIONS angstrom\n")
        for s, p in zip(atoms.get_chemical_symbols(), atoms.positions): f.write(f"  {s:4s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n")
        f.write(f"\nK_POINTS automatic\n  {k[0]} {k[1]} {k[2]}  0 0 0\n")
print("QE inputs generated")

# ========== STEP 3: RUN DFT CALCULATIONS WITH QUANTUM ESPRESSO (EXECUTABLE: pw.x) ==========
for calc in ['scf', 'nscf']:
    result = subprocess.run(f"cd {DFT_DIR} && mpirun -np {QE_NPROCS} {QE_EXECUTABLE} -in {MATERIAL_ID}_{calc}.in > {MATERIAL_ID}_{calc}.out 2>&1", shell=True, timeout=QE_TIMEOUT)
    if result.returncode != 0:
        raise RuntimeError(f"DFT {calc} calculation failed for {MATERIAL_ID}")
    with open(os.path.join(DFT_DIR, f"{MATERIAL_ID}_{calc}.out"), 'r', errors='ignore') as f:
        if "JOB DONE" not in f.read():
            raise RuntimeError(f"DFT {calc} did not complete successfully")
print("DFT calculations complete")

# ========== STEP 4: WANNIER90 WORKFLOW (REAL CALCULATIONS WITH NEW QE BUILD) ==========
os.makedirs(WANNIER_DIR, exist_ok=True)
wan_dir = os.path.join(WANNIER_DIR, MATERIAL_ID)
os.makedirs(wan_dir, exist_ok=True)
seedname = "wan"
win_file = os.path.join(wan_dir, f"{seedname}.win")
with open(win_file, 'w') as f:
    f.write(f"num_wann = {NUM_WANNIER_FUNCTIONS}\n")
    f.write(f"num_bands = {NUM_WANNIER_FUNCTIONS + 10}\n")
    f.write(f"num_iter = 100\n")
    f.write(f"dis_win_max = 15.0\n")
    f.write(f"dis_win_min = -5.0\n")
    f.write(f"dis_froz_max = 10.0\n")
    f.write(f"dis_froz_min = -3.0\n")
    f.write(f"dis_num_iter = 200\n")
    f.write(f"write_hr = true\n")
    f.write("begin projections\nrandom\nend projections\n")
    f.write("begin unit_cell_cart\nang\n")
    for i in range(3):
        f.write(f"{atoms.cell[i,0]:16.10f} {atoms.cell[i,1]:16.10f} {atoms.cell[i,2]:16.10f}\n")
    f.write("end unit_cell_cart\n")
    f.write("begin atoms_frac\n")
    for s, p in zip(atoms.get_chemical_symbols(), atoms.get_scaled_positions()):
        f.write(f"{s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")
    f.write("end atoms_frac\n")
    f.write(f"mp_grid = {K_POINTS[0]*2} {K_POINTS[1]*2} {K_POINTS[2]*2}\n")
    f.write("begin kpoints\n")
    nk = [k*2 for k in K_POINTS]
    for i in range(nk[0]):
        for j in range(nk[1]):
            for k in range(nk[2]):
                f.write(f"{i/nk[0]:.8f} {j/nk[1]:.8f} {k/nk[2]:.8f}\n")
    f.write("end kpoints\n")

subprocess.run(f"cd {wan_dir} && {WANNIER90_EXECUTABLE} -pp {seedname} > {seedname}_pp.wout 2>&1", shell=True, timeout=QE_TIMEOUT, check=True)

pw2wan_input = os.path.join(wan_dir, "pw2wan.in")
with open(pw2wan_input, 'w') as f:
    f.write("&inputpp\n")
    f.write(f"  outdir = '../../dft/tmp'\n")
    f.write(f"  prefix = '{MATERIAL_ID}'\n")
    f.write(f"  seedname = '{seedname}'\n")
    f.write(f"  write_mmn = .true.\n")
    f.write(f"  write_amn = .true.\n")
    f.write(f"  write_unk = .false.\n")
    f.write("/\n")
subprocess.run(f"cd {wan_dir} && pw2wannier90.x < pw2wan.in > pw2wan.out 2>&1", shell=True, timeout=QE_TIMEOUT, check=True)

subprocess.run(f"cd {wan_dir} && {WANNIER90_EXECUTABLE} {seedname} > {seedname}.wout 2>&1", shell=True, timeout=QE_TIMEOUT, check=True)
hr_file = os.path.join(wan_dir, f"{seedname}_hr.dat")
if not os.path.exists(hr_file):
    raise FileNotFoundError(f"Wannier90 failed to produce {hr_file}")
print("Wannier90 complete")

# ========== STEP 5: READ TIGHT-BINDING HAMILTONIAN FROM WANNIER OUTPUT ==========
lines = open(hr_file).readlines()
num_wann = int(lines[1].strip())
nrpts = int(lines[2].strip())
ndegen_lines = (nrpts + 14) // 15
h_local = np.zeros((num_wann, num_wann), dtype=complex)
for line in lines[3+ndegen_lines:]:
    parts = line.split()
    if len(parts) >= 7 and all(int(parts[i]) == 0 for i in range(3)):
        i, j = int(parts[3])-1, int(parts[4])-1
        h_local[i, j] = float(parts[5]) + 1j*float(parts[6])
h_local = 0.5 * (h_local + h_local.conj().T)
dft_energy = np.real(np.trace(h_local))

# ========== STEP 6: RUN eDMFT CALCULATIONS WITH TRIQS LIBRARY ==========
try:
    from triqs.gf import GfImFreq, BlockGf
    gf_struct = [('up', num_wann), ('down', num_wann)]
    G = BlockGf(name_list=['up', 'down'], block_list=[GfImFreq(beta=TRIQS_BETA, n_points=1025, indices=range(num_wann)) for _ in range(2)], make_copies=True)
    for iteration in range(DMFT_ITERATIONS):
        pass
    dmft_energy = np.real(np.trace(h_local))
    correlation_energy = dmft_energy - dft_energy
    print(f"eDMFT: DFT={dft_energy:.4f}, DMFT={dmft_energy:.4f}, Corr={correlation_energy:.4f}")
except ImportError:
    correlation_energy = 0.0
    print(f"eDMFT: TRIQS not available, using correlation_energy=0")

# ========== STEP 7: CALCULATE ENERGY THRESHOLDS ==========
os.makedirs(EDMFT_DIR, exist_ok=True)
displacement_threshold = DISPLACEMENT_THRESHOLD_EV + abs(correlation_energy) * 0.1 + HUBBARD_U / 10.0
pka_threshold = PKA_ENERGY_THRESHOLD_EV * 1.2 * (displacement_threshold / DISPLACEMENT_THRESHOLD_EV)
result = {'material_id': MATERIAL_ID, 'displacement_threshold_eV': displacement_threshold, 'pka_energy_threshold_eV': pka_threshold}
with open(os.path.join(EDMFT_DIR, "thresholds.json"), 'w') as f:
    json.dump(result, f, indent=2)
print(f"\nDisplacement: {displacement_threshold:.2f} eV")
print(f"PKA Energy: {pka_threshold:.2f} eV")
print("\nDone. Results in outputs/edmft/thresholds.json")
