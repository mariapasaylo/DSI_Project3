#!/usr/bin/env python3
import os
import sys
import subprocess


def check_dependencies():
    packages = {'jarvis': 'jarvis-tools', 'ase': 'ase', 'numpy': 'numpy'}
    missing = []
    
    for module, package in packages.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(package)
    
    return missing


def test_step1():
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'pipeline'))
    
    try:
        import step1_material_selection
        import config
        
        config.N_MATERIALS = 1
        results = step1_material_selection.main()
        
        return results and len(results) > 0
    except Exception as e:
        print(f"Error: {e}")
        return False


def main():
    print(f"{'='*60}")
    print("Pipeline Test")
    print(f"{'='*60}\n")
    
    missing = check_dependencies()
    if missing:
        print(f"Missing: {', '.join(missing)}")
        return False
    
    print("✓ Dependencies OK")
    
    print("\nTesting material selection...")
    if test_step1():
        print("✓ Material selection works\n")
        print("Ready to run: cd pipeline && python run_pipeline.py")
        return True
    else:
        print("✗ Material selection failed")
        return False


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
