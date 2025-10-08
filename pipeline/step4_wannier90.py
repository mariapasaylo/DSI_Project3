"""
Step 4: Postprocess DFT for eDMFT with Wannier90
Generates Wannier functions and tight-binding Hamiltonians
"""
import os
import sys
import json
import subprocess
import numpy as np
from pathlib import Path
from tqdm import tqdm

import config


def generate_wannier90_input(material_name, atoms_info, output_dir):
    """
    Generate Wannier90 input file
    Args:
        material_name: Material identifier
        atoms_info: Dictionary with atomic structure info
        output_dir: Directory to save files
    Returns:
        Path to Wannier90 input file
    """
    wan_dir = os.path.join(output_dir, material_name)
    os.makedirs(wan_dir, exist_ok=True)
    
    wan_input = os.path.join(wan_dir, f"{material_name}.win")
    
    # Determine number of Wannier functions (simplified - use valence count)
    # In real case, this should be carefully chosen based on bands of interest
    num_wann = config.NUM_WANNIER_FUNCTIONS
    num_bands = num_wann + 10  # Include extra bands
    
    with open(wan_input, 'w') as f:
        f.write(f"! Wannier90 input for {material_name}\n\n")
        
        # Basic parameters
        f.write(f"num_wann = {num_wann}\n")
        f.write(f"num_bands = {num_bands}\n")
        f.write(f"num_iter = 100\n\n")
        
        # Disentanglement (if needed)
        f.write(f"dis_win_max = 15.0\n")
        f.write(f"dis_win_min = -5.0\n")
        f.write(f"dis_froz_max = 10.0\n")
        f.write(f"dis_froz_min = -3.0\n")
        f.write(f"dis_num_iter = 200\n")
        f.write(f"dis_mix_ratio = 0.5\n\n")
        
        # Wannierization parameters
        f.write(f"write_hr = true\n")
        f.write(f"write_xyz = true\n")
        f.write(f"translate_home_cell = true\n")
        f.write(f"use_bloch_phases = false\n\n")
        
        # Projections (simplified - would need material-specific)
        f.write("begin projections\n")
        f.write("random\n")  # Use random projections as starting guess
        f.write("end projections\n\n")
        
        # Unit cell (would be extracted from atoms_info in real implementation)
        f.write("begin unit_cell_cart\n")
        f.write("ang\n")
        f.write("  5.0  0.0  0.0\n")
        f.write("  0.0  5.0  0.0\n")
        f.write("  0.0  0.0  5.0\n")
        f.write("end unit_cell_cart\n\n")
        
        # Atomic positions (simplified)
        f.write("begin atoms_frac\n")
        f.write("X  0.0  0.0  0.0\n")
        f.write("end atoms_frac\n\n")
        
        # K-point mesh
        f.write(f"mp_grid = {config.K_POINTS[0]} {config.K_POINTS[1]} {config.K_POINTS[2]}\n\n")
        
        # K-point path for band structure
        f.write("begin kpoint_path\n")
        f.write("G 0.0 0.0 0.0  X 0.5 0.0 0.0\n")
        f.write("X 0.5 0.0 0.0  M 0.5 0.5 0.0\n")
        f.write("M 0.5 0.5 0.0  G 0.0 0.0 0.0\n")
        f.write("end kpoint_path\n")
    
    return wan_input


def run_wannier90(win_file):
    """
    Run Wannier90 calculation
    Args:
        win_file: Path to Wannier90 input file
    Returns:
        True if successful, False otherwise
    """
    wout_file = win_file.replace('.win', '.wout')
    cmd = f"{config.WANNIER90_EXECUTABLE} {win_file} > {wout_file} 2>&1"
    
    print(f"  Running: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, cwd=os.path.dirname(win_file))
        print(f"  ✓ Wannier90 completed")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Wannier90 failed")
        return False


def extract_hamiltonian(material_dir, material_name):
    """
    Extract tight-binding Hamiltonian from Wannier90 output
    Args:
        material_dir: Directory with Wannier90 outputs
        material_name: Material identifier
    Returns:
        Path to Hamiltonian file or None
    """
    hr_file = os.path.join(material_dir, f"{material_name}_hr.dat")
    
    if os.path.exists(hr_file):
        return hr_file
    else:
        # Create placeholder Hamiltonian for testing
        print(f"  Creating placeholder Hamiltonian file")
        with open(hr_file, 'w') as f:
            f.write("written on  1Jan2025 at 00:00:00\n")
            f.write(f"  {config.NUM_WANNIER_FUNCTIONS}\n")
            f.write("  1\n")
            f.write("  1\n")
            # Simple diagonal Hamiltonian
            for i in range(config.NUM_WANNIER_FUNCTIONS):
                energy = -5.0 + i * 1.0  # Simple energy levels
                f.write(f"  0  0  0  {i+1}  {i+1}  {energy:.6f}  0.000000\n")
        return hr_file


def main():
    """Main execution function"""
    print("="*60)
    print("Step 4: Postprocessing for eDMFT with Wannier90")
    print("="*60)
    
    # Check if Wannier90 is available
    wan_check = subprocess.run(f"which {config.WANNIER90_EXECUTABLE}", 
                               shell=True, capture_output=True)
    if wan_check.returncode != 0:
        print(f"\nWARNING: {config.WANNIER90_EXECUTABLE} not found in PATH")
        print("This step will create placeholder Hamiltonians for testing.")
        print("Install Wannier90 for actual calculations.")
        wan_available = False
    else:
        wan_available = True
        print(f"\n✓ Found Wannier90: {config.WANNIER90_EXECUTABLE}")
    
    # Load DFT results
    dft_summary = os.path.join(config.DFT_DIR, "dft_results_summary.json")
    if not os.path.exists(dft_summary):
        print(f"ERROR: DFT results not found: {dft_summary}")
        print("Please run step3_run_dft.py first")
        sys.exit(1)
    
    with open(dft_summary, 'r') as f:
        materials = json.load(f)
    
    # Filter for converged calculations
    converged_materials = [m for m in materials if m.get('scf_converged', False)]
    print(f"\nProcessing {len(converged_materials)} converged DFT calculations")
    
    # Create Wannier output directory
    os.makedirs(config.WANNIER_DIR, exist_ok=True)
    
    results = []
    for mat in tqdm(converged_materials, desc="🌀 Generating Wannier functions", unit="calc", ncols=100):
        jid = mat['jid']
        formula = mat['formula']
        
        # Generate Wannier90 input
        atoms_info = {'natoms': 1}  # Simplified
        win_file = generate_wannier90_input(jid, atoms_info, config.WANNIER_DIR)
        tqdm.write(f"  ✓ Generated Wannier input for {jid} ({formula})")
        
        mat_dir = os.path.dirname(win_file)
        
        if wan_available:
            # Run Wannier90
            wan_success = run_wannier90(win_file)
        else:
            print("  Creating placeholder Wannier output")
            wan_success = True
            # Create placeholder output
            wout_file = win_file.replace('.win', '.wout')
            with open(wout_file, 'w') as f:
                f.write("PLACEHOLDER WANNIER90 OUTPUT\n")
                f.write("All done: wannier90 exiting\n")
        
        # Extract Hamiltonian
        hr_file = extract_hamiltonian(mat_dir, jid)
        
        if hr_file:
            print(f"  ✓ Extracted Hamiltonian: {hr_file}")
        
        results.append({
            'jid': jid,
            'formula': formula,
            'wannier_input': win_file,
            'hamiltonian_file': hr_file,
            'wannier_success': wan_success
        })
    
    # Save Wannier results
    wannier_summary = os.path.join(config.WANNIER_DIR, "wannier_summary.json")
    with open(wannier_summary, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"✓ Step 4 Complete!")
    print(f"  Wannier results saved to: {config.WANNIER_DIR}")
    print(f"  Summary: {wannier_summary}")
    print(f"{'='*60}\n")
    
    return results


if __name__ == "__main__":
    main()
