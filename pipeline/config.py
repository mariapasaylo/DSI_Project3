import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPELINE_DIR = os.path.join(BASE_DIR, "pipeline")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
MATERIALS_DIR = os.path.join(OUTPUT_DIR, "materials")
DFT_DIR = os.path.join(OUTPUT_DIR, "dft")
WANNIER_DIR = os.path.join(OUTPUT_DIR, "wannier")
EDMFT_DIR = os.path.join(OUTPUT_DIR, "edmft")

N_MATERIALS = 1
RANDOM_SEED = None
JARVIS_DATABASE = "dft_3d"

QE_PSEUDOPOTENTIALS_DIR = "/usr/share/espresso/pseudo"
QE_EXECUTABLE = "pw.x"
QE_NPROCS = 8
QE_TIMEOUT = 600
QE_MAX_ITERATIONS = 50

ECUTWFC = 30.0
ECUTRHO = 120.0
K_POINTS = [2, 2, 2]
OCCUPATIONS = "smearing"
SMEARING = "gaussian"
DEGAUSS = 0.02

WANNIER90_EXECUTABLE = "wannier90.x"
NUM_WANNIER_FUNCTIONS = 10

TRIQS_BETA = 40.0
TRIQS_N_ITER = 20
HUBBARD_U = 4.0
CHEMICAL_POTENTIAL = 0.0

ASK_CLEAR_OUTPUTS = True
