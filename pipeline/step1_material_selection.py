import os
import sys
import json
import random
from tqdm import tqdm

from jarvis.db.figshare import data as jarvis_data
from jarvis.core.atoms import Atoms as JarvisAtoms
from ase.io import write as ase_write
from ase import Atoms as AseAtoms

import config


def jarvis_to_ase(jarvis_atoms):
    return AseAtoms(
        symbols=jarvis_atoms.elements,
        positions=jarvis_atoms.cart_coords,
        cell=jarvis_atoms.lattice_mat,
        pbc=True
    )


def select_random_materials(n=1, database="dft_3d"):
    if config.RANDOM_SEED is not None:
        random.seed(config.RANDOM_SEED)
    
    print(f"Loading JARVIS {database}...")
    all_data = jarvis_data(database)
    valid = [e for e in all_data if 'atoms' in e and 'jid' in e]
    print(f"Found {len(valid)} materials")
    
    return random.sample(valid, min(n, len(valid)))


def save_structure(material_data, output_dir):
    jid = material_data['jid']
    formula = material_data.get('formula', 'unknown')
    
    mat_dir = os.path.join(output_dir, jid)
    os.makedirs(mat_dir, exist_ok=True)
    
    jarvis_atoms = JarvisAtoms.from_dict(material_data['atoms'])
    ase_atoms = jarvis_to_ase(jarvis_atoms)
    
    files = {
        'cif': os.path.join(mat_dir, f"{jid}.cif"),
        'poscar': os.path.join(mat_dir, f"{jid}.vasp"),
        'xyz': os.path.join(mat_dir, f"{jid}.xyz"),
        'json': os.path.join(mat_dir, f"{jid}_metadata.json")
    }
    
    ase_write(files['cif'], ase_atoms)
    ase_write(files['poscar'], ase_atoms, format='vasp')
    ase_write(files['xyz'], ase_atoms)
    
    metadata = {
        'jid': jid,
        'formula': formula,
        'natoms': len(ase_atoms),
        'elements': list(set(ase_atoms.get_chemical_symbols())),
        'lattice_parameters': ase_atoms.cell.cellpar().tolist(),
        'jarvis_data': {k: v for k, v in material_data.items() if k != 'atoms'}
    }
    
    with open(files['json'], 'w') as f:
        json.dump(metadata, f, indent=2)
    
    return files, metadata


def main():
    os.makedirs(config.MATERIALS_DIR, exist_ok=True)
    
    materials = select_random_materials(
        n=config.N_MATERIALS,
        database=config.JARVIS_DATABASE
    )
    
    results = []
    for mat in tqdm(materials, desc="Processing", unit="mat", ncols=80):
        jid = mat['jid']
        formula = mat.get('formula', 'unknown')
        files, metadata = save_structure(mat, config.MATERIALS_DIR)
        
        tqdm.write(f"✓ {jid} ({formula}) - {metadata['natoms']} atoms")
        
        results.append({
            'jid': jid,
            'formula': formula,
            'files': files,
            'metadata': metadata
        })
    
    summary_file = os.path.join(config.MATERIALS_DIR, "materials_summary.json")
    with open(summary_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    return results


if __name__ == "__main__":
    main()
