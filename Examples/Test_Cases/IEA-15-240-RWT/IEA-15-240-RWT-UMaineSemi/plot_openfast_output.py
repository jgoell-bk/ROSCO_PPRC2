#!/usr/bin/env python3
"""
Plot OpenFAST output variables over time
Usage: python plot_openfast_output.py <variable_name> [--tmax TIME]
Example: python plot_openfast_output.py PtfmPitch --tmax 80
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def read_openfast_output(outfile):
    """Read OpenFAST .out file into a pandas DataFrame"""
    # Read all lines
    with open(outfile, 'r') as f:
        lines = f.readlines()
    
    # Find the header line (contains "Time" and other column names)
    header_idx = None
    for i, line in enumerate(lines):
        if 'Time' in line and 's)' not in line.lower():  # Column name line (not units)
            header_idx = i
            break
    
    if header_idx is None:
        raise ValueError("Could not find column headers")
    
    # Parse column names and handle duplicates
    col_names = lines[header_idx].split()
    
    # Make column names unique by appending numbers to duplicates
    seen = {}
    unique_names = []
    for name in col_names:
        if name in seen:
            seen[name] += 1
            unique_names.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            unique_names.append(name)
    
    # Skip header and units line, read data
    data = pd.read_csv(outfile, sep=r'\s+', skiprows=header_idx+2, 
                       names=unique_names, low_memory=False)
    
    return data

def plot_variable(outfile, variable_name, tmax=None):
    """
    Plot a specific variable from OpenFAST output file
    
    Parameters:
    -----------
    outfile : str
        Path to the .out file
    variable_name : str
        Name of the variable to plot (e.g., 'PtfmPitch', 'GenTq', 'GenPwr')
    tmax : float, optional
        Maximum time to plot (seconds). If None, plots all data.
    """
    # Read the output file
    print(f"Reading {outfile}...")
    data = read_openfast_output(outfile)
    
    # Filter by time if tmax is specified
    if tmax is not None:
        data = data[data['Time'] <= tmax]
        print(f"Filtering data to t <= {tmax} seconds")
    
    # Check if variable exists
    if variable_name not in data.columns:
        print(f"\nError: Variable '{variable_name}' not found in output file.")
        print(f"\nAvailable variables:")
        for i, col in enumerate(data.columns, 1):
            print(f"  {col}")
            if i % 3 == 0:  # Print 3 columns per line
                print()
        return
    
    # Get time and variable data
    time = data['Time']
    
    # Find the variable (try with and without units in brackets)
    var_col = None
    for col in data.columns:
        if variable_name in col or col.startswith(variable_name):
            var_col = col
            break
    
    if var_col is None:
        print(f"\nError: Variable '{variable_name}' not found in output file.")
        print(f"\nAvailable variables:")
        for i, col in enumerate(data.columns, 1):
            print(f"  {col}")
            if i % 3 == 0:  # Print 3 columns per line
                print()
        return
    
    variable = data[var_col]
    
    # Get units from column name (format is typically "VarName_(units)")
    units = var_col.split('(')[-1].rstrip(')') if '(' in var_col else ''
    var_name_clean = var_col.split('_(')[0] if '_(' in var_col else var_col
    
    # Create plot
    plt.figure(figsize=(12, 6))
    plt.plot(time, variable, linewidth=1.5)
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel(f'{var_name_clean} [{units}]' if units else var_name_clean, fontsize=12)
    plt.title(f'{var_name_clean} vs Time', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save and show
    output_png = outfile.replace('.out', f'_{var_name_clean}.png')
    plt.savefig(output_png, dpi=150)
    print(f"\nPlot saved to: {output_png}")
    plt.show()


def plot_multiple_variables(outfile, variable_names, tmax=None):
    """
    Plot multiple variables on subplots
    
    Parameters:
    -----------
    outfile : str
        Path to the .out file
    variable_names : list
        List of variable names to plot
    tmax : float, optional
        Maximum time to plot (seconds). If None, plots all data.
    """
    # Read the output file
    print(f"Reading {outfile}...")
    data = read_openfast_output(outfile)
    
    # Filter by time if tmax is specified
    if tmax is not None:
        data = data[data['Time'] <= tmax]
        print(f"Filtering data to t <= {tmax} seconds")
    
    # Filter valid variables
    valid_vars = []
    for var_name in variable_names:
        var_col = None
        for col in data.columns:
            if var_name in col or col.startswith(var_name):
                var_col = col
                break
        if var_col:
            valid_vars.append(var_col)
        else:
            print(f"Warning: Variable '{var_name}' not found")
    
    if not valid_vars:
        print("No valid variables to plot!")
        return
    
    # Get time data
    time = data['Time']
    
    # Create subplots
    n_vars = len(valid_vars)
    fig, axes = plt.subplots(n_vars, 1, figsize=(12, 4*n_vars))
    if n_vars == 1:
        axes = [axes]
    
    for ax, var_col in zip(axes, valid_vars):
        variable = data[var_col]
        
        # Get units and clean name
        units = var_col.split('(')[-1].rstrip(')') if '(' in var_col else ''
        var_name_clean = var_col.split('_(')[0] if '_(' in var_col else var_col
        
        ax.plot(time, variable, linewidth=1.5)
        ax.set_xlabel('Time [s]', fontsize=11)
        ax.set_ylabel(f'{var_name_clean} [{units}]' if units else var_name_clean, fontsize=11)
        ax.set_title(f'{var_name_clean}', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save and show
    output_png = outfile.replace('.out', '_multi.png')
    plt.savefig(output_png, dpi=150)
    print(f"\nPlot saved to: {output_png}")
    plt.show()


if __name__ == "__main__":
    # Default output file
    outfile = "IEA-15-240-RWT-UMaineSemi.out"
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Plot OpenFAST output variables',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot_openfast_output.py PtfmPitch
  python plot_openfast_output.py PtfmPitch --tmax 80
  python plot_openfast_output.py BldPitch1 BldPitch2 BldPitch3 --tmax 60
  python plot_openfast_output.py PtfmPitch GenTq GenPwr

Common PPPR-relevant variables:
  - PtfmPitch      : Platform pitch angle
  - GenTq          : Generator torque
  - GenPwr         : Generator power
  - GenSpeed       : Generator speed
  - BldPitch1/2/3  : Blade pitch angles
  - Wind1VelX      : Wind speed
        """
    )
    parser.add_argument('variables', nargs='+', help='Variable name(s) to plot')
    parser.add_argument('--tmax', type=float, default=None, 
                        help='Maximum time to plot (seconds). Default: plot all data')
    
    args = parser.parse_args()
    
    # Plot single or multiple variables
    if len(args.variables) == 1:
        plot_variable(outfile, args.variables[0], tmax=args.tmax)
    else:
        plot_multiple_variables(outfile, args.variables, tmax=args.tmax)