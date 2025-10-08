"""
Configuration file for DFT → eDMFT Pipeline
Centralized settings for all workflow components
"""
import os

# Base directories
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPELINE_DIR = os.path.join(BASE_DIR, "pipeline")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
MATERIALS_DIR = os.path.join(OUTPUT_DIR, "materials")
DFT_DIR = os.path.join(OUTPUT_DIR, "dft")
WANNIER_DIR = os.path.join(OUTPUT_DIR, "wannier")
EDMFT_DIR = os.path.join(OUTPUT_DIR, "edmft")

# Material selection settings
N_MATERIALS = 1  # Number of random materials to select (reduced for testing)
JARVIS_DATABASE = "dft_3d"  # JARVIS database to query

# Quantum ESPRESSO settings
QE_PSEUDOPOTENTIALS_DIR = "/usr/share/espresso/pseudo"  # System QE pseudos
QE_EXECUTABLE = "pw.x"  # Adjust if QE is not in PATH
QE_NPROCS = 8  # Number of processors for parallel QE runs (USE ALL CORES!)
QE_TIMEOUT = 600  # Max time per calculation (seconds) - 10 minutes
QE_MAX_ITERATIONS = 50  # Max SCF iterations before giving up

# DFT calculation parameters
ECUTWFC = 30.0  # Wavefunction cutoff (Ry) - REDUCED for speed
ECUTRHO = 120.0  # Charge density cutoff (Ry) - REDUCED for speed
K_POINTS = [2, 2, 2]  # K-point grid - REDUCED for speed (was [4,4,4])
OCCUPATIONS = "smearing"
SMEARING = "gaussian"
DEGAUSS = 0.02

# Wannier90 settings
WANNIER90_EXECUTABLE = "wannier90.x"
NUM_WANNIER_FUNCTIONS = 10  # Will be adjusted per material

# eDMFT/TRIQS settings
TRIQS_BETA = 40.0  # Inverse temperature (1/eV)
TRIQS_N_ITER = 20  # Number of DMFT iterations
HUBBARD_U = 4.0  # Hubbard U parameter (eV)
CHEMICAL_POTENTIAL = 0.0  # Initial guess

# Output settings
SAVE_INTERMEDIATE_FILES = True
VERBOSE = True
