#!/usr/bin/env python3
"""
Quick test script to validate pipeline installation and basic functionality
"""
import os
import sys

def check_python_packages():
    """Check if required Python packages are installed"""
    print("Checking Python packages...")
    packages = {
        'jarvis': 'jarvis-tools',
        'ase': 'ase',
        'pymatgen': 'pymatgen',
        'numpy': 'numpy',
        'scipy': 'scipy',
    }
    
    missing = []
    for module, package in packages.items():
        try:
            __import__(module)
            print(f"  ✓ {package}")
        except ImportError:
            print(f"  ✗ {package} (missing)")
            missing.append(package)
    
    # Check TRIQS separately (optional)
    try:
        import triqs
        print(f"  ✓ triqs")
    except ImportError:
        print(f"  ⚠ triqs (optional - will use placeholder mode)")
    
    return missing


def check_external_tools():
    """Check if external tools are available"""
    print("\nChecking external tools...")
    import subprocess
    
    tools = {
        'pw.x': 'Quantum ESPRESSO',
        'wannier90.x': 'Wannier90',
    }
    
    missing = []
    for cmd, name in tools.items():
        result = subprocess.run(f"which {cmd}", shell=True, 
                               capture_output=True, text=True)
        if result.returncode == 0:
            print(f"  ✓ {name} ({cmd})")
        else:
            print(f"  ⚠ {name} (not in PATH - will use placeholder mode)")
            missing.append(name)
    
    return missing


def test_imports():
    """Test importing pipeline modules"""
    print("\nTesting pipeline imports...")
    
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'pipeline'))
    
    try:
        import config
        print(f"  ✓ config")
    except Exception as e:
        print(f"  ✗ config: {e}")
        return False
    
    modules = [
        'step1_material_selection',
        'step2_convert_to_qe',
        'step3_run_dft',
        'step4_wannier90',
        'step5_edmft',
        'run_pipeline'
    ]
    
    for mod in modules:
        try:
            __import__(mod)
            print(f"  ✓ {mod}")
        except Exception as e:
            print(f"  ✗ {mod}: {e}")
            return False
    
    return True


def test_step1():
    """Test Step 1 - Material selection"""
    print("\nTesting Step 1 (Material Selection)...")
    
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'pipeline'))
        import step1_material_selection
        import config
        
        # Override to select just 1 material for quick test
        original_n = config.N_MATERIALS
        config.N_MATERIALS = 1
        
        print("  Running material selection (1 material)...")
        results = step1_material_selection.main()
        
        config.N_MATERIALS = original_n
        
        if results and len(results) > 0:
            print(f"  ✓ Step 1 successful - selected {len(results)} material(s)")
            return True
        else:
            print(f"  ✗ Step 1 failed - no results")
            return False
            
    except Exception as e:
        print(f"  ✗ Step 1 error: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests"""
    print("="*60)
    print("  DFT → eDMFT Pipeline Test Suite")
    print("="*60 + "\n")
    
    # Check dependencies
    missing_packages = check_python_packages()
    missing_tools = check_external_tools()
    
    # Test imports
    import_ok = test_imports()
    
    if not import_ok:
        print("\n" + "="*60)
        print("⚠ Pipeline imports failed - check installation")
        print("="*60)
        return False
    
    # Test Step 1 (quick functional test)
    if missing_packages:
        print("\n" + "="*60)
        print("⚠ Missing required packages:")
        for pkg in missing_packages:
            print(f"    pip install {pkg}")
        print("="*60)
        return False
    
    step1_ok = test_step1()
    
    print("\n" + "="*60)
    print("  Test Results Summary")
    print("="*60)
    
    if missing_packages:
        print(f"✗ Missing packages: {', '.join(missing_packages)}")
        print(f"  Run: pip install {' '.join(missing_packages)}")
    else:
        print("✓ All required Python packages installed")
    
    if missing_tools:
        print(f"⚠ External tools not in PATH: {', '.join(missing_tools)}")
        print("  Pipeline will use placeholder mode for these steps")
    else:
        print("✓ All external tools available")
    
    if import_ok:
        print("✓ All pipeline modules import successfully")
    else:
        print("✗ Pipeline module import failed")
    
    if step1_ok:
        print("✓ Step 1 functional test passed")
    else:
        print("✗ Step 1 functional test failed")
    
    print("\n" + "="*60)
    
    if import_ok and not missing_packages:
        print("✓ READY TO RUN PIPELINE")
        print("\nRun: cd pipeline && python run_pipeline.py")
    else:
        print("⚠ INSTALL MISSING DEPENDENCIES")
    
    print("="*60 + "\n")
    
    return import_ok and not missing_packages


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
