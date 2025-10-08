#!/usr/bin/env python3
import os
import sys
import shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import config
import step1_material_selection
import step2_convert_to_qe
import step3_run_dft
import step4_wannier90
import step5_edmft


def clear_outputs():
    if not config.ASK_CLEAR_OUTPUTS:
        return
    
    if not os.path.exists(config.OUTPUT_DIR):
        return
    
    response = input(f"\nClear previous outputs in {config.OUTPUT_DIR}? (y/N): ").strip().lower()
    if response == 'y':
        shutil.rmtree(config.OUTPUT_DIR)
        print(f"✓ Cleared {config.OUTPUT_DIR}")


def run_full_pipeline(start_step=1, end_step=5):
    clear_outputs()
    
    print(f"\n{'='*60}")
    print("DFT → eDMFT Pipeline")
    print(f"{'='*60}\n")
    
    steps = {
        1: ("Material Selection", step1_material_selection.main),
        2: ("QE Input Generation", step2_convert_to_qe.main),
        3: ("DFT Calculation", step3_run_dft.main),
        4: ("Wannier90", step4_wannier90.main),
        5: ("eDMFT (TRIQS)", step5_edmft.main),
    }
    
    for step_num in range(start_step, end_step + 1):
        step_name, step_func = steps[step_num]
        print(f"\n[Step {step_num}] {step_name}")
        print("-" * 60)
        
        try:
            step_func()
            print(f"✓ Step {step_num} complete\n")
        except Exception as e:
            print(f"✗ Step {step_num} failed: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    
    print(f"\n{'='*60}")
    print("Pipeline Complete")
    print(f"{'='*60}")
    print(f"Output: {os.path.join(config.EDMFT_DIR, 'niel_calculator_input.json')}\n")


if __name__ == "__main__":
    run_full_pipeline()
