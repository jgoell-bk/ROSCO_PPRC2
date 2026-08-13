#!/usr/bin/env python3
"""
Plot OpenFAST output variables over time
Usage: python plot_openfast_output.py <variable_name> [--tmax TIME]
Example: python plot_openfast_output.py PtfmPitch --tmax 80
"""

import sys
import numpy as np
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


def _find_col(data, name):
    """Resolve a base variable name to its actual column (handles _N duplicates and unit suffixes)."""
    for col in data.columns:
        if col == name or col.startswith(name + '_(') or col == name + '_1':
            return col
    for col in data.columns:
        if name in col:
            return col
    return None


def plot_pppr_diagnostic(outfile, tmax=None,
                         pppr_amp_phi_deg=5.0,
                         pppr_freq_hz=0.03391,
                         pppr_activation_time=30.0,
                         pc_maxpit_deg=90.0,
                         vs_maxtq_knm=19790.0,
                         rated_genspeed_rpm=7.55):
    """
    PPPR diagnostic dashboard for the PI: 2x3 panels showing the failure cascade.

    Panels:
      (1) BlPitch1 with PC_MaxPit saturation line
      (2) PtfmPitch with sinusoidal phi reference overlay
      (3) GenSpeed with rated speed line
      (4) GenTq with VS_MaxTq saturation line
      (5) TwrBsMyt (tower base fore-aft moment)
      (6) RotThrust
    """
    print(f"Reading {outfile}...")
    data = read_openfast_output(outfile)

    if tmax is not None:
        data = data[data['Time'] <= tmax]
        print(f"Filtering data to t <= {tmax} seconds")

    t = data['Time'].values

    fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharex=True)

    # (1) BlPitch1 with saturation line
    ax = axes[0, 0]
    col = _find_col(data, 'BldPitch1')
    if col is not None:
        ax.plot(t, data[col], linewidth=1.2, color='C0')
    ax.axhline(pc_maxpit_deg, color='red', linestyle='--', linewidth=1,
               label=f'PC_MaxPit = {pc_maxpit_deg:.0f}°')
    ax.set_ylabel('BlPitch1 [deg]')
    ax.set_title('Blade pitch (saturation = full feather)')
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)

    # (2) PtfmPitch with sinusoidal reference overlay
    ax = axes[0, 1]
    col = _find_col(data, 'PtfmPitch')
    if col is not None:
        ax.plot(t, data[col], linewidth=1.2, color='C0', label='PtfmPitch (measured)')
    phi_ref = np.where(
        t >= pppr_activation_time,
        pppr_amp_phi_deg * np.sin(2 * np.pi * pppr_freq_hz * t),
        np.nan,
    )
    ax.plot(t, phi_ref, linewidth=1.0, color='red', linestyle='--',
            label=f'phi_ref ({pppr_amp_phi_deg:.1f}° @ {pppr_freq_hz:.4f} Hz)')
    ax.axvline(pppr_activation_time, color='gray', linestyle=':', linewidth=0.8,
               label=f'PPPR on (t={pppr_activation_time:.0f}s)')
    ax.set_ylabel('PtfmPitch [deg]')
    ax.set_title('Platform pitch vs reference')
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)

    # (3) GenSpeed with rated speed line
    ax = axes[0, 2]
    col = _find_col(data, 'GenSpeed')
    if col is not None:
        ax.plot(t, data[col], linewidth=1.2, color='C0')
    ax.axhline(rated_genspeed_rpm, color='green', linestyle='--', linewidth=1,
               label=f'Rated = {rated_genspeed_rpm:.2f} rpm')
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_ylabel('GenSpeed [rpm]')
    ax.set_title('Generator speed (negative = rotor reversal)')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    # (4) GenTq with VS_MaxTq saturation line
    ax = axes[1, 0]
    col = _find_col(data, 'GenTq')
    if col is not None:
        ax.plot(t, data[col], linewidth=1.2, color='C0')
    ax.axhline(vs_maxtq_knm, color='red', linestyle='--', linewidth=1,
               label=f'VS_MaxTq = {vs_maxtq_knm:.0f} kN·m')
    ax.set_ylabel('GenTq [kN·m]')
    ax.set_xlabel('Time [s]')
    ax.set_title('Generator torque (saturation = controller windup)')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    # (5) TwrBsMyt - tower base fore-aft moment
    ax = axes[1, 1]
    col = _find_col(data, 'TwrBsMyt')
    if col is not None:
        ax.plot(t, data[col], linewidth=1.2, color='C0')
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_ylabel('TwrBsMyt [kN·m]')
    ax.set_xlabel('Time [s]')
    ax.set_title('Tower base fore-aft moment (structural impact)')
    ax.grid(True, alpha=0.3)

    # (6) RotThrust
    ax = axes[1, 2]
    col = _find_col(data, 'RotThrust')
    if col is not None:
        ax.plot(t, data[col], linewidth=1.2, color='C0')
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_ylabel('RotThrust [kN]')
    ax.set_xlabel('Time [s]')
    ax.set_title('Rotor thrust (collapse on feather)')
    ax.grid(True, alpha=0.3)

    fig.suptitle('PPPR Diagnostic: failure cascade overview',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    output_png = outfile.replace('.out', '_pppr_diagnostic.png')
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
  python plot_openfast_output.py --pppr                       # 2x3 PPPR diagnostic dashboard
  python plot_openfast_output.py --pppr --tmax 200            # diagnostic, first 200 s
  python plot_openfast_output.py --pppr --pppr-amp-phi 5.0 --pppr-freq 0.03391

Common PPPR-relevant variables:
  - PtfmPitch      : Platform pitch angle
  - GenTq          : Generator torque
  - GenPwr         : Generator power
  - GenSpeed       : Generator speed
  - BldPitch1/2/3  : Blade pitch angles
  - Wind1VelX      : Wind speed
        """
    )
    parser.add_argument('variables', nargs='*', help='Variable name(s) to plot (omit when using --pppr)')
    parser.add_argument('--tmax', type=float, default=None,
                        help='Maximum time to plot (seconds). Default: plot all data')
    parser.add_argument('--pppr', action='store_true',
                        help='Render the 2x3 PPPR diagnostic dashboard')
    parser.add_argument('--pppr-amp-phi', type=float, default=5.0,
                        help='Phi reference amplitude in degrees (default 5.0)')
    parser.add_argument('--pppr-freq', type=float, default=0.03391,
                        help='PPPR reference frequency in Hz (default 0.03391 = 0.213 rad/s)')
    parser.add_argument('--pppr-activation', type=float, default=30.0,
                        help='PPPR activation time in seconds (default 30.0)')
    parser.add_argument('--pc-maxpit', type=float, default=90.0,
                        help='Pitch saturation upper limit in deg (default 90)')
    parser.add_argument('--vs-maxtq', type=float, default=19790.0,
                        help='Torque saturation limit in kN·m (default 19790)')
    parser.add_argument('--rated-genspeed', type=float, default=7.55,
                        help='Rated generator speed in rpm (default 7.55)')

    args = parser.parse_args()

    if args.pppr:
        plot_pppr_diagnostic(
            outfile,
            tmax=args.tmax,
            pppr_amp_phi_deg=args.pppr_amp_phi,
            pppr_freq_hz=args.pppr_freq,
            pppr_activation_time=args.pppr_activation,
            pc_maxpit_deg=args.pc_maxpit,
            vs_maxtq_knm=args.vs_maxtq,
            rated_genspeed_rpm=args.rated_genspeed,
        )
    elif len(args.variables) == 1:
        plot_variable(outfile, args.variables[0], tmax=args.tmax)
    elif len(args.variables) > 1:
        plot_multiple_variables(outfile, args.variables, tmax=args.tmax)
    else:
        parser.error('Provide variable names, or use --pppr for the diagnostic dashboard.')