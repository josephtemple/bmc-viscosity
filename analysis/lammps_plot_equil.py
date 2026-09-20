"""
Plot equilibration-important variables (potential energy, temperature, pressure) as a function of time, from csv
Note: this NEEDS to be run AFTER lammpslog_to_csv.py, as that creates the csv this file reads from

Usage:
    python lammps_plot_equil.py run_name

Expects a thermo_style like:
    thermo_style custom step pe press ke temp lx ly lz pxx pyy pzz spcpu density
"""
import os
import sys
import csv
import matplotlib
import pandas as pd
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot_thermo(df, out_png):
    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "stix"
    })

    steps = df['Step']

    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)

    axes[0].plot(steps, df['PotEng'], color="tab:blue")
    axes[0].set_ylabel("Potential Energy (kcal/mol)")
    axes[0].set_title("Thermo diagnostics")

    axes[1].plot(steps, df['Temp'], color="tab:red")
    axes[1].set_ylabel("Temperature (K)")

    axes[2].plot(steps, df['Press'], color="tab:green", label="Press (avg)")
    
    axes[2].plot(steps, df['Pxx'], alpha=0.4, label="Pxx")
    axes[2].plot(steps, df['Pyy'], alpha=0.4, label="Pyy")
    axes[2].plot(steps, df['Pzz'], alpha=0.4, label="Pzz")

    axes[2].set_ylabel("Pressure (atm)")
    axes[2].set_xlabel("Timestep (fs)")
    axes[2].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print(f"saved {out_png}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python run_equil.py run_name")
        sys.exit(1)

    file_name = sys.argv[1] + ".csv"
    png_name = sys.argv[1] + ".png"

    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv_dir = os.path.join(script_dir, "equil_results")
    csv_file = os.path.join(csv_dir, file_name)

    png_file = os.path.join(csv_dir, png_name)

    df = pd.read_csv(csv_file, header=0)

    plot_thermo(df, png_file)