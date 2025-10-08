import os
import sys
import json
import glob
from tqdm import tqdm
from ase.io import read as ase_read

import config


PSEUDO_MAP = {
    'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
    'Li': 'Li.pbe-sl-rrkjus_psl.1.0.0.UPF',
    'C': 'C.pbe-n-rrkjus_psl.1.0.0.UPF',
    'N': 'N.pbe-n-rrkjus_psl.1.0.0.UPF',
    'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
    'F': 'F.pbe-n-rrkjus_psl.1.0.0.UPF',
    'Na': 'Na.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Mg': 'Mg.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Al': 'Al.pbe-n-rrkjus_psl.1.0.0.UPF',
    'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
    'P': 'P.pbe-n-rrkjus_psl.1.0.0.UPF',
    'S': 'S.pbe-n-rrkjus_psl.1.0.0.UPF',
    'Cl': 'Cl.pbe-n-rrkjus_psl.1.0.0.UPF',
    'K': 'K.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Ca': 'Ca.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Ti': 'Ti.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'V': 'V.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Cr': 'Cr.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Mn': 'Mn.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Co': 'Co.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Ni': 'Ni.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Cu': 'Cu.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Zn': 'Zn.pbe-dnl-rrkjus_psl.1.0.0.UPF',
    'Ga': 'Ga.pbe-dnl-rrkjus_psl.1.0.0.UPF',
    'Ge': 'Ge.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'As': 'As.pbe-n-rrkjus_psl.1.0.0.UPF',
    'Se': 'Se.pbe-n-rrkjus_psl.1.0.0.UPF',
    'Br': 'Br.pbe-n-rrkjus_psl.1.0.0.UPF',
}


def get_pseudopotentials(elements):
    pseudos = {}
    pseudo_dir = config.QE_PSEUDOPOTENTIALS_DIR
    
    for elem in elements:
        if elem in PSEUDO_MAP:
            pseudo_file = PSEUDO_MAP[elem]
            if os.path.exists(os.path.join(pseudo_dir, pseudo_file)):
                pseudos[elem] = pseudo_file
                continue
        
        patterns = [
            f"{elem}.pbe-*.UPF",
            f"{elem}.paw*.upf",
            f"{elem}_pbe*.UPF",
            f"{elem.lower()}_pbe*.UPF",
            f"{elem}_ONCV*.upf",
            f"{elem}*.oncvpsp.upf",
            f"{elem}.*.UPF",
            f"{elem}*.upf",
        ]
        
        found = False
        for pattern in patterns:
            matches = glob.glob(os.path.join(pseudo_dir, pattern))
            if matches:
                matches_sorted = sorted(matches, key=lambda x: (
                    'paw' not in x.lower(),
                    'oncv' not in x.lower(),
                    'pbe' not in x.lower(),
                ))
                pseudos[elem] = os.path.basename(matches_sorted[0])
                found = True
                break
        
        if not found:
            pseudos[elem] = f"{elem}.pbe-n-rrkjus_psl.1.0.0.UPF"
    
    return pseudos


def generate_qe_scf_input(atoms, material_name, output_dir):
    elements = list(set(atoms.get_chemical_symbols()))
    pseudos = get_pseudopotentials(elements)
    
    qe_dir = os.path.join(output_dir, material_name)
    os.makedirs(qe_dir, exist_ok=True)
    scf_input = os.path.join(qe_dir, f"{material_name}_scf.in")
    
    with open(scf_input, 'w') as f:
        f.write("&CONTROL\n")
        f.write(f"  calculation = 'scf'\n")
        f.write(f"  prefix = '{material_name}'\n")
        f.write(f"  outdir = './tmp'\n")
        f.write(f"  pseudo_dir = '{config.QE_PSEUDOPOTENTIALS_DIR}'\n")
        f.write(f"  verbosity = 'high'\n")
        f.write(f"  wf_collect = .false.\n")
        f.write(f"/\n\n")
        
        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0\n")
        f.write(f"  nat = {len(atoms)}\n")
        f.write(f"  ntyp = {len(elements)}\n")
        f.write(f"  ecutwfc = {config.ECUTWFC}\n")
        f.write(f"  ecutrho = {config.ECUTRHO}\n")
        f.write(f"  occupations = '{config.OCCUPATIONS}'\n")
        f.write(f"  smearing = '{config.SMEARING}'\n")
        f.write(f"  degauss = {config.DEGAUSS}\n")
        f.write(f"/\n\n")
        
        # Electrons section  
        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = 1.0d-6\n")  # Relaxed for speed
        f.write(f"  mixing_beta = 0.7\n")
        f.write(f"  electron_maxstep = {getattr(config, 'QE_MAX_ITERATIONS', 50)}\n")
        f.write(f"  diagonalization = 'david'\n")  # Davidson for parallel efficiency
        f.write(f"/\n\n")
        
        # Atomic species
        f.write("ATOMIC_SPECIES\n")
        # Rough mass estimates (would be better to use actual values)
        atomic_masses = {elem: atoms.get_masses()[atoms.get_chemical_symbols().index(elem)] 
                        for elem in elements}
        for elem in elements:
            f.write(f"  {elem}  {atomic_masses[elem]:.4f}  {pseudos[elem]}\n")
        f.write("\n")
        
        # Cell parameters
        f.write("CELL_PARAMETERS angstrom\n")
        cell = atoms.get_cell()
        for i in range(3):
            f.write(f"  {cell[i, 0]:16.10f} {cell[i, 1]:16.10f} {cell[i, 2]:16.10f}\n")
        f.write("\n")
        
        # Atomic positions
        f.write("ATOMIC_POSITIONS angstrom\n")
        for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            f.write(f"  {symbol:4s} {pos[0]:16.10f} {pos[1]:16.10f} {pos[2]:16.10f}\n")
        f.write("\n")
        
        # K-points
        f.write("K_POINTS automatic\n")
        f.write(f"  {config.K_POINTS[0]} {config.K_POINTS[1]} {config.K_POINTS[2]}  0 0 0\n")
    
    return scf_input


def generate_qe_nscf_input(atoms, material_name, output_dir):
    elements = list(set(atoms.get_chemical_symbols()))
    pseudos = get_pseudopotentials(elements)
    
    qe_dir = os.path.join(output_dir, material_name)
    nscf_input = os.path.join(qe_dir, f"{material_name}_nscf.in")
    k_dense = [k * 2 for k in config.K_POINTS]
    
    with open(nscf_input, 'w') as f:
        f.write("&CONTROL\n")
        f.write(f"  calculation = 'nscf'\n")
        f.write(f"  prefix = '{material_name}'\n")
        f.write(f"  outdir = './tmp'\n")
        f.write(f"  pseudo_dir = '{config.QE_PSEUDOPOTENTIALS_DIR}'\n")
        f.write(f"  verbosity = 'high'\n")
        f.write(f"/\n\n")
        
        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0\n")
        f.write(f"  nat = {len(atoms)}\n")
        f.write(f"  ntyp = {len(elements)}\n")
        f.write(f"  ecutwfc = {config.ECUTWFC}\n")
        f.write(f"  ecutrho = {config.ECUTRHO}\n")
        f.write(f"  occupations = '{config.OCCUPATIONS}'\n")
        f.write(f"  smearing = '{config.SMEARING}'\n")
        f.write(f"  degauss = {config.DEGAUSS}\n")
        f.write(f"  nosym = .true.\n")
        f.write(f"  noinv = .true.\n")
        f.write(f"/\n\n")
        
        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = 1.0d-10\n")
        f.write(f"  diago_full_acc = .true.\n")
        f.write(f"/\n\n")
        
        f.write("ATOMIC_SPECIES\n")
        atomic_masses = {elem: atoms.get_masses()[atoms.get_chemical_symbols().index(elem)] 
                        for elem in elements}
        for elem in elements:
            f.write(f"  {elem}  {atomic_masses[elem]:.4f}  {pseudos[elem]}\n")
        f.write("\n")
        
        f.write("CELL_PARAMETERS angstrom\n")
        cell = atoms.get_cell()
        for i in range(3):
            f.write(f"  {cell[i, 0]:16.10f} {cell[i, 1]:16.10f} {cell[i, 2]:16.10f}\n")
        f.write("\n")
        
        f.write("ATOMIC_POSITIONS angstrom\n")
        for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            f.write(f"  {symbol:4s} {pos[0]:16.10f} {pos[1]:16.10f} {pos[2]:16.10f}\n")
        f.write("\n")
        
        f.write("K_POINTS automatic\n")
        f.write(f"  {k_dense[0]} {k_dense[1]} {k_dense[2]}  0 0 0\n")
    
    return nscf_input


def main():
    summary_file = os.path.join(config.MATERIALS_DIR, "materials_summary.json")
    if not os.path.exists(summary_file):
        print(f"ERROR: Materials summary not found")
        sys.exit(1)
    
    with open(summary_file, 'r') as f:
        materials = json.load(f)
    
    os.makedirs(config.DFT_DIR, exist_ok=True)
    
    results = []
    for mat in tqdm(materials, desc="Generating", unit="mat", ncols=80):
        jid = mat['jid']
        poscar_file = mat['files']['poscar']
        atoms = ase_read(poscar_file)
        
        scf_input = generate_qe_scf_input(atoms, jid, config.DFT_DIR)
        nscf_input = generate_qe_nscf_input(atoms, jid, config.DFT_DIR)
        
        tqdm.write(f"✓ {jid}")
        
        results.append({
            'jid': jid,
            'formula': mat['formula'],
            'scf_input': scf_input,
            'nscf_input': nscf_input
        })
    
    conversion_summary = os.path.join(config.DFT_DIR, "qe_inputs_summary.json")
    with open(conversion_summary, 'w') as f:
        json.dump(results, f, indent=2)
    
    return results


if __name__ == "__main__":
    main()
