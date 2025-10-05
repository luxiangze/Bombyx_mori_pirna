#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Add scripts directory to Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the fixed plot_pirna_trim function
from scripts.pirna_analysis.visualization import plot_pirna_trim

def main():
    """Test the fixed plot_pirna_trim function"""
    # Set input and output directories
    input_dir = "/home/gyk/project/ld_pirna/results/pirna_analysis_20250405_164942"
    output_dir = "/home/gyk/project/ld_pirna/results/test_fixed_plot_pirna_trim"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read pirna_trim data for all samples
    all_pirna_trim_data = {}
    sample_files = [f for f in os.listdir(input_dir) if f.endswith("_pirna_trim.csv")]
    
    print(f"Found the following sample files:")
    for sample_file in sample_files:
        sample_name = sample_file.replace("_pirna_trim.csv", "")
        print(f"  - {sample_file}")
        
        # Read CSV file
        file_path = os.path.join(input_dir, sample_file)
        df = pd.read_csv(file_path)
        
        # Print basic information
        print(f"    Rows: {len(df)}")
        print(f"    Columns: {', '.join(df.columns)}")
        
        # Store data
        all_pirna_trim_data[sample_name] = df
    
    print(f"\nTotal samples found: {len(all_pirna_trim_data)}")
    
    # Call the fixed plot_pirna_trim function
    print(f"\nStart plotting piRNA trimming comparison...")
    plot_pirna_trim(all_pirna_trim_data, f"{output_dir}/")
    
    print(f"\nPlotting completed. Results saved to: {output_dir}")

if __name__ == "__main__":
    main()
