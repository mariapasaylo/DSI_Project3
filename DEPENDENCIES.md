# Complete Dependency List - DSI Project 3
## DFT+eDMFT Pipeline for Space Photovoltaics

Last Updated: October 16, 2025

---

## System Environment

### Operating System
- **OS**: Ubuntu 25.04 (Plucky Puffin)
- **Architecture**: ARM64/aarch64 (Apple Silicon compatible)
- **Kernel**: Linux
- **Link**: https://ubuntu.com/

### Package Manager
- **apt**: System package manager for Ubuntu
- **conda**: 25.9.0 (Miniconda3)
- **pip**: Python package installer

---

## Python Environment

### Python Interpreter
- **Version**: Python 3.9.23 (DSI conda environment)
- **Base Python**: Python 3.13.5 (Miniconda3 base)
- **Conda Environment**: `DSI` (at `/home/vm/miniconda3/envs/DSI/`)
- **Links**: 
  - https://www.python.org/
  - https://docs.conda.io/

### Python Packages (DSI environment)

#### Core Scientific Computing
1. **numpy**: 2.0.2
   - Numerical arrays and linear algebra
   - https://numpy.org/

2. **scipy**: 1.13.1
   - Scientific computing algorithms
   - https://scipy.org/

3. **pandas**: 2.3.3
   - Data analysis and manipulation
   - https://pandas.pydata.org/

#### Visualization
4. **matplotlib**: 3.9.4
   - 2D plotting library
   - https://matplotlib.org/

5. **seaborn**: 0.13.2
   - Statistical data visualization
   - https://seaborn.pydata.org/

6. **nglview**: 4.0
   - 3D molecular structure visualization (WebGL)
   - https://github.com/nglviewer/nglview

#### Atomic Simulation
7. **ase**: 3.26.0 (Atomic Simulation Environment)
   - Atomic structure manipulation and I/O
   - Used for material structure handling
   - https://wiki.fysik.dtu.dk/ase/

#### Materials Database
8. **jarvis-tools**: 2025.5.30
   - NIST JARVIS materials database access
   - Access to 3D crystal structures and properties
   - https://jarvis.nist.gov/
   - https://github.com/usnistgov/jarvis

#### Interactive Computing
9. **IPython**: 8.18.1
   - Enhanced Python shell
   - https://ipython.org/

10. **Jupyter**: 1.1.1
    - Jupyter Notebook system
    - https://jupyter.org/

11. **JupyterLab**: 4.4.9
    - Next-gen Jupyter interface
    - https://jupyterlab.readthedocs.io/

---

## Scientific Computing Software

### Quantum ESPRESSO (DFT Calculations)

**Version**: 6.7-3build2  
**Source**: Ubuntu 25.04 ARM64 repository  
**Installation**: `apt install quantum-espresso quantum-espresso-data-sssp`  
**Location**: `/usr/bin/`

**Key Executables**:
- `pw.x` - Plane-Wave self-consistent field (SCF) and non-SCF calculations
- `pw2wannier90.x` - Interface to Wannier90 for Wannier function generation
- `ph.x` - Phonon calculations
- `pp.x` - Post-processing tools
- `dos.x` - Density of states
- Plus 70+ additional utilities

**Linked Libraries** (ARM64 native):
- ScaLAPACK 2.2 (distributed linear algebra)
- LAPACK 3 (linear algebra)
- OpenBLAS 3 (optimized BLAS)
- FFTW3 (Fast Fourier Transform)
- Open MPI 4.1.1 (MPI implementation)
- gfortran 5 (Fortran runtime)

**Pseudopotentials**:
- Package: `quantum-espresso-data-sssp` 1.3.0-3
- Library: Standard Solid-State Pseudopotentials (SSSP)
- Location: `/usr/share/espresso/pseudo/`
- Format: UPF (Unified Pseudopotential Format)

**Configuration**:
- Parallel execution: 8 MPI processes
- Timeout: 600 seconds
- Cutoff energy: 30 eV (wavefunction), 120 eV (charge density)
- K-points: 2×2×2 automatic grid (SCF), 4×4×4 (NSCF)

**Links**:
- https://www.quantum-espresso.org/
- https://github.com/QEF/q-e
- https://www.materialscloud.org/discover/sssp (pseudopotentials)

---

### Wannier90 (Wannier Function Construction)

**Version**: 3.1.0+ds-10  
**Source**: Ubuntu 25.04 ARM64 repository  
**Installation**: `apt install wannier90`  
**Location**: `/usr/bin/wannier90.x` (854 KB binary)

**Purpose**: Constructs maximally-localized Wannier functions from Bloch states

**Build Configuration**:
- Compiler: gfortran
- Optimization: -O3
- Build reference: `/home/vm/wannier90_build/wannier90-3.1.0/`

**Linked Libraries**:
- LAPACK (linear algebra)
- BLAS (basic linear algebra)
- gfortran runtime

**Configuration**:
- Number of Wannier functions: 10
- Projections: Random initial guess
- Iterations: 100
- Output: Hamiltonian in real space (wan_hr.dat)

**Links**:
- https://wannier.org/
- https://github.com/wannier-developers/wannier90
- https://wannier.readthedocs.io/

---

### TRIQS (Quantum Many-Body Physics)

**Version**: 3.x (built from source)  
**Source**: Built from source for ARM64  
**Build Location**: `/home/vm/triqs_build/triqs/`  
**Install Location**: `/home/vm/miniconda3/envs/DSI/triqs_install/`

**Full Name**: Toolbox for Research on Interacting Quantum Systems

**Purpose**: 
- Many-body quantum physics framework
- Green's function calculations
- Dynamical mean-field theory (DMFT)
- Quantum Monte Carlo solvers

**Components Installed**:
- Core TRIQS library (`libtriqs.so`)
- Python bindings (in DSI site-packages)
- C++ headers
- CMake configuration files

**Build Configuration**:
- CMake build system
- Python support: Enabled (Python 3.9.23)
- NumPy integration: Enabled
- Install prefix: `/home/vm/miniconda3/envs/DSI/triqs_install/`
- Compiler: gcc/g++
- Build type: Release

**Linked Libraries**:
- HDF5 (data storage)
- BLAS/LAPACK (linear algebra)
- MPI (parallel computing)
- FFTW (Fourier transforms)
- Python 3.9 (bindings)

**Python Modules**:
- `triqs` - Core framework
- `triqs.gf` - Green's functions (GfImFreq, BlockGf)
- `triqs.operators` - Second quantization operators
- `triqs.version` - Version information

**Links**:
- https://triqs.github.io/
- https://github.com/TRIQS/triqs
- https://triqs.github.io/triqs/latest/
- https://triqs.github.io/triqs/latest/documentation/manual/

---

### TRIQS CTHYB (QMC Impurity Solver)

**Version**: 3.x (built from source)  
**Source**: Built from source for ARM64 (requires TRIQS)  
**Install Location**: `/home/vm/miniconda3/envs/DSI/triqs_install/`

**Full Name**: Continuous-Time Hybridization-expansion Quantum Monte Carlo solver

**Purpose**:
- Solves quantum impurity problems in DMFT
- Uses continuous-time QMC with hybridization expansion
- Captures quantum fluctuations and correlations beyond mean-field

**Components**:
- `triqs_cthyb` Python module
- `libtriqs_cthyb_c.a` - C++ library
- CMake configuration files

**Build Configuration**:
- Depends on: TRIQS 3.x
- Compiler: gcc/g++
- Build type: Release
- MPI support: Yes (via TRIQS)

**QMC Parameters** (typical):
- Monte Carlo cycles: 10,000
- Cycle length: 100
- Warmup cycles: 5,000
- Measurements: Green's function, density, self-energy

**Python API**:
- `triqs_cthyb.Solver` - Main solver class
- Integration with TRIQS Green's functions
- Operator algebra from `triqs.operators`

**Links**:
- https://triqs.github.io/cthyb/
- https://github.com/TRIQS/cthyb
- https://triqs.github.io/cthyb/latest/
- Tutorial: https://github.com/TRIQS/tutorials/tree/main/ModelDMFT

---

### Open MPI (Message Passing Interface)

**Version**: 4.1.1 (HYDRA build, March 2023)  
**Source**: Ubuntu 25.04 repository  
**Executables**: `mpirun`, `mpiexec`  
**Location**: `/usr/bin/`

**Purpose**: Parallel computing framework for distributed memory systems

**Used By**:
- Quantum ESPRESSO (8 parallel processes)
- Wannier90 (optional parallel postprocessing)
- TRIQS (parallel QMC calculations)

**Components**:
- `libmpi.so.40` - MPI library
- `libucp.so.0` - Unified communication protocol
- `libucs.so.0` - Unified communication services
- `libfabric.so.1` - Fabric abstraction
- `libevent` - Event notification
- `libhwloc.so.15` - Hardware locality

**Configuration**:
- Default: 8 processes (`mpirun -np 8`)
- Network: UCX transport layer
- Hardware: ARM64 optimized

**Link**: https://www.open-mpi.org/

---

## System Libraries (ARM64)

### Linear Algebra
- **OpenBLAS**: Optimized BLAS implementation for ARM64
- **LAPACK**: Linear Algebra PACKage
- **ScaLAPACK**: Scalable LAPACK for distributed computing

### Numerical Libraries
- **FFTW3**: Fastest Fourier Transform in the West
- **HDF5**: Hierarchical Data Format (used by TRIQS)

### Compiler Runtimes
- **gfortran 5**: GNU Fortran runtime library
- **libgcc**: GCC runtime library
- **libstdc++**: C++ standard library

### System Libraries
- **glibc**: GNU C Library
- **libm**: Math library
- **libpthread**: POSIX threads

---

## Build Tools & Compilers

### Compilers
- **gcc/g++**: GNU C/C++ compiler (for TRIQS)
- **gfortran**: GNU Fortran compiler (for QE, Wannier90)

### Build Systems
- **CMake**: Cross-platform build system (for TRIQS)
- **Make**: GNU Make (for Wannier90)
- **autotools**: For QE configuration

### Version Control
- **git**: Used to clone TRIQS source

---

## Jupyter Ecosystem

### Core Packages
- **jupyter**: 1.1.1 - Metapackage
- **jupyter-client**: 8.6.3 - Jupyter protocol client
- **jupyter-core**: 5.8.1 - Core Jupyter functionality
- **jupyter-server**: 2.17.0 - Jupyter server

### JupyterLab
- **jupyterlab**: 4.4.9 - Lab interface
- **jupyterlab-server**: 2.27.3 - Lab server
- **jupyterlab-widgets**: 3.0.15 - Widget extensions
- **jupyterlab-pygments**: 0.3.0 - Syntax highlighting

### Notebook
- **jupyter-console**: 6.6.3 - Terminal console
- **notebook**: Package for classic notebook interface

### Extensions
- **jupyter-events**: 0.12.0 - Event system
- **jupyter-lsp**: 2.3.0 - Language Server Protocol
- **jupyter-server-terminals**: 0.5.3 - Terminal support

---

## Standard Library Modules (Python)

Used throughout the notebook:
- `os` - Operating system interface
- `sys` - System-specific parameters
- `re` - Regular expressions
- `subprocess` - Process management
- `tempfile` - Temporary file creation
- `random` - Random number generation
- `time` - Time measurement
- `importlib` - Dynamic imports
- `shutil` - High-level file operations

---

## Pipeline Configuration

### DFT Parameters
- **Plane-wave cutoff**: 30.0 eV (wavefunction)
- **Charge density cutoff**: 120.0 eV (4× wavefunction)
- **K-point grid**: 2×2×2 (SCF), 4×4×4 (NSCF)
- **Occupation**: Smearing (Gaussian, σ=0.02)
- **Convergence**: 1.0×10⁻⁸ (electrons)

### Wannier90 Parameters
- **Number of functions**: 10
- **Bands**: 30 (from NSCF)
- **Iterations**: 100
- **Projections**: Random initial guess

### DMFT Parameters
- **Hubbard U**: 4.0 eV
- **Inverse temperature (β)**: 40.0 eV⁻¹
- **DMFT iterations**: 20
- **Orbitals**: 1 (single-orbital model)
- **Chemical potential**: Calculated from DFT

### QMC Parameters
- **Monte Carlo cycles**: 10,000
- **Cycle length**: 100 steps
- **Warmup cycles**: 5,000
- **Measurements**: Every cycle

---

## File Structure

### Working Directory
- **WORK_DIR**: `/home/vm/DSI_Project3/calculations/`

### Calculation Directories
- `dft_{MATERIAL_ID}/` - DFT calculations
  - `{MATERIAL_ID}_scf.in` - SCF input
  - `{MATERIAL_ID}_nscf.in` - NSCF input
  - `tmp/` - QE scratch directory
  
- `wannier_{MATERIAL_ID}/` - Wannier90 calculations
  - `wan.win` - Wannier90 input
  - `wan_hr.dat` - Hamiltonian in real space
  - `wan.amn`, `wan.mmn`, `wan.eig` - Overlap matrices

---

## Data Sources

### JARVIS Database
- **Database**: `dft_3d` (3D materials)
- **Source**: NIST JARVIS (Joint Automated Repository for Various Integrated Simulations)
- **Materials**: 40,000+ DFT-calculated materials
- **Properties**: Formation energy, band gap, elastic constants, etc.
- **Link**: https://jarvis.nist.gov/

### Test Material
- **ID**: JVASP-25075
- **Formula**: Te₂ (Tellurium)
- **Source**: JARVIS DFT database

---

## Build Instructions Summary

### TRIQS Build (from source)
```bash
cd /home/vm/triqs_build/triqs
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/home/vm/miniconda3/envs/DSI/triqs_install \
         -DPYTHON_EXECUTABLE=/home/vm/miniconda3/envs/DSI/bin/python
make -j8
make install
```

### TRIQS CTHYB Build (from source)
```bash
# After TRIQS is installed
cd /home/vm/cthyb_build
cmake . -DCMAKE_INSTALL_PREFIX=/home/vm/miniconda3/envs/DSI/triqs_install
make -j8
make install
```

### System Packages (Ubuntu ARM64)
```bash
sudo apt install quantum-espresso quantum-espresso-data-sssp wannier90
```

---

## Known Limitations

1. **ARM64 Architecture**: 
   - TRIQS and CTHYB must be built from source (no pre-built ARM64 binaries)
   - QE and Wannier90 available as Ubuntu ARM64 packages

2. **Python Environment**:
   - TRIQS requires Python 3.9 (DSI environment)
   - Notebook runs in DSI conda environment

3. **Computational Resources**:
   - DFT calculations: ~30-60 seconds per material (8 cores)
   - Wannier90: ~10-30 seconds
   - DMFT QMC: ~5-10 minutes per iteration (20 iterations)

---

## References

### Software Documentation
- Quantum ESPRESSO: Giannozzi et al., J. Phys.: Condens. Matter 21, 395502 (2009)
- Wannier90: Pizzi et al., J. Phys.: Condens. Matter 32, 165902 (2020)
- TRIQS: Parcollet et al., Comput. Phys. Commun. 196, 398-415 (2015)
- TRIQS CTHYB: Seth et al., Comput. Phys. Commun. 200, 274-284 (2016)

### Scientific Methods
- MLWF: Marzari & Vanderbilt, Phys. Rev. B 56, 12847 (1997)
- DMFT: Georges et al., Rev. Mod. Phys. 68, 13 (1996)
- SSSP: Prandini et al., npj Comput. Mater. 4, 72 (2018)

---

**Total Dependencies**: 11 Python packages + 4 major software suites + system libraries + build tools
