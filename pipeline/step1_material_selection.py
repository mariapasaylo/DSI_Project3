"""
Step 1: Material Selection & Structure Input from JARVIS Dataset
Fetches crystal structures from JARVIS and saves them in multiple formats
"""
import os
import sys
import json
import random
from pathlib import Path
from tqdm import tqdm

try:
    from jarvis.db.figshare import data as jarvis_data
    from jarvis.core.atoms import Atoms as JarvisAtoms
except ImportError:
    print("ERROR: jarvis-tools not installed. Run: pip install jarvis-tools")
    sys.exit(1)

try:
    from ase.io import write as ase_write
    from ase import Atoms as AseAtoms
except ImportError:
    print("ERROR: ase not installed. Run: pip install ase")
    sys.exit(1)

import config


def jarvis_to_ase(jarvis_atoms):
    """Convert JARVIS Atoms to ASE Atoms object"""
    lattice = jarvis_atoms.lattice_mat
    positions = jarvis_atoms.cart_coords
    symbols = jarvis_atoms.elements
    
    ase_atoms = AseAtoms(
        symbols=symbols,
        positions=positions,
        cell=lattice,
        pbc=True
    )
    return ase_atoms


def select_random_materials(n=2, database="dft_3d"):
    """
    Select random materials from JARVIS database
    Args:
        n: Number of materials to select
        database: JARVIS database name
    Returns:
        List of material dictionaries
    """
    print(f"Loading JARVIS {database} database...")
    
    # Load JARVIS data
    all_data = jarvis_data(database)
    
    # Filter for materials with complete data
    valid_materials = []
    for entry in all_data:
        if 'atoms' in entry and 'jid' in entry:
            valid_materials.append(entry)
    
    print(f"Found {len(valid_materials)} valid materials in database")
    
    # Select random materials
    selected = random.sample(valid_materials, min(n, len(valid_materials)))
    
    return selected


def save_structure(material_data, output_dir):
    """
    Save material structure in multiple formats
    Args:
        material_data: JARVIS material dictionary
        output_dir: Directory to save files
    Returns:
        Dictionary with file paths
    """
    jid = material_data['jid']
    formula = material_data.get('formula', 'unknown')
    
    # Create material-specific directory
    mat_dir = os.path.join(output_dir, jid)
    os.makedirs(mat_dir, exist_ok=True)
    
    # Convert JARVIS atoms to ASE
    jarvis_atoms = JarvisAtoms.from_dict(material_data['atoms'])
    ase_atoms = jarvis_to_ase(jarvis_atoms)
    
    # Save in multiple formats
    files = {}
    files['cif'] = os.path.join(mat_dir, f"{jid}.cif")
    files['poscar'] = os.path.join(mat_dir, f"{jid}.vasp")
    files['xyz'] = os.path.join(mat_dir, f"{jid}.xyz")
    files['json'] = os.path.join(mat_dir, f"{jid}_metadata.json")
    
    ase_write(files['cif'], ase_atoms)
    ase_write(files['poscar'], ase_atoms, format='vasp')
    ase_write(files['xyz'], ase_atoms)
    
    # Save metadata
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
    """Main execution function"""
    print("="*60)
    print("Step 1: Material Selection from JARVIS Dataset")
    print("="*60)
    
    # Create output directory
    os.makedirs(config.MATERIALS_DIR, exist_ok=True)
    
    # Select materials
    materials = select_random_materials(
        n=config.N_MATERIALS,
        database=config.JARVIS_DATABASE
    )
    
    print(f"\nSelected {len(materials)} materials:")
    
    # Process and save each material
    results = []
    for mat in tqdm(materials, desc="📦 Processing materials", unit="material", ncols=100):
        jid = mat['jid']
        formula = mat.get('formula', 'unknown')
        
        files, metadata = save_structure(mat, config.MATERIALS_DIR)
        
        tqdm.write(f"  ✓ {jid} ({formula}) - {', '.join(metadata['elements'])} [{metadata['natoms']} atoms]")
        
        results.append({
            'jid': jid,
            'formula': formula,
            'files': files,
            'metadata': metadata
        })
    
    # Save summary
    summary_file = os.path.join(config.MATERIALS_DIR, "materials_summary.json")
    with open(summary_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"✓ Step 1 Complete!")
    print(f"  Materials saved to: {config.MATERIALS_DIR}")
    print(f"  Summary: {summary_file}")
    print(f"{'='*60}\n")
    
    return results


if __name__ == "__main__":
    main()
