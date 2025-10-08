#!/usr/bin/env python3
"""
Master Pipeline Runner
Executes complete DFT → eDMFT workflow autonomously
"""
import os
import sys
import argparse
from pathlib import Path

# Add pipeline directory to path
pipeline_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, pipeline_dir)

import config
import step1_material_selection
import step2_convert_to_qe
import step3_run_dft
import step4_wannier90
import step5_edmft


def run_full_pipeline(start_step=1, end_step=5):
    """
    Run complete pipeline from material selection to eDMFT
    Args:
        start_step: Starting step number (1-5)
        end_step: Ending step number (1-5)
    """
    print("\n" + "="*70)
    print("  DFT → eDMFT PIPELINE")
    print("  Autonomous Workflow for Correlated Materials Simulation")
    print("="*70 + "\n")
    
    steps = {
        1: ("Material Selection from JARVIS", step1_material_selection.main),
        2: ("Convert to Quantum ESPRESSO Input", step2_convert_to_qe.main),
        3: ("Run DFT Simulation", step3_run_dft.main),
        4: ("Wannier90 Postprocessing", step4_wannier90.main),
        5: ("eDMFT Calculation with TRIQS", step5_edmft.main),
    }
    
    for step_num in range(start_step, end_step + 1):
        if step_num not in steps:
            print(f"ERROR: Invalid step number {step_num}")
            continue
        
        step_name, step_func = steps[step_num]
        
        print(f"\n{'='*70}")
        print(f"  EXECUTING STEP {step_num}: {step_name}")
        print(f"{'='*70}\n")
        
        try:
            result = step_func()
            print(f"\n✓ Step {step_num} completed successfully\n")
        except Exception as e:
            print(f"\n✗ ERROR in Step {step_num}: {e}")
            print(f"Pipeline halted at step {step_num}\n")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    
    print("\n" + "="*70)
    print("  PIPELINE COMPLETE!")
    print("="*70)
    print("\nOutput directories:")
    print(f"  Materials:     {config.MATERIALS_DIR}")
    print(f"  DFT Results:   {config.DFT_DIR}")
    print(f"  Wannier Data:  {config.WANNIER_DIR}")
    print(f"  eDMFT Results: {config.EDMFT_DIR}")
    print(f"\nNIEL Calculator Input: {os.path.join(config.EDMFT_DIR, 'niel_calculator_input.json')}")
    print("\nReady for NIEL calculator integration!")
    print("="*70 + "\n")


def main():
    """Main entry point with command-line interface"""
    parser = argparse.ArgumentParser(
        description='DFT → eDMFT Pipeline for Correlated Materials',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  python run_pipeline.py
  
  # Run specific steps
  python run_pipeline.py --start 1 --end 3
  
  # Run single step
  python run_pipeline.py --start 2 --end 2
  
Pipeline Steps:
  1. Material Selection from JARVIS dataset
  2. Convert structures to Quantum ESPRESSO input
  3. Run DFT simulation with QE
  4. Generate Wannier functions with Wannier90
  5. Run eDMFT calculations with TRIQS
        """
    )
    
    parser.add_argument(
        '--start', 
        type=int, 
        default=1, 
        choices=[1, 2, 3, 4, 5],
        help='Starting step (default: 1)'
    )
    
    parser.add_argument(
        '--end', 
        type=int, 
        default=5, 
        choices=[1, 2, 3, 4, 5],
        help='Ending step (default: 5)'
    )
    
    parser.add_argument(
        '--config',
        type=str,
        help='Path to custom config file (optional)'
    )
    
    args = parser.parse_args()
    
    if args.start > args.end:
        print("ERROR: Start step must be <= end step")
        sys.exit(1)
    
    # Load custom config if provided
    if args.config:
        print(f"Loading custom config from {args.config}")
        # Would implement config loading here
    
    # Run pipeline
    run_full_pipeline(args.start, args.end)


if __name__ == "__main__":
    main()
