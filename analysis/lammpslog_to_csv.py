"""
Parse thermo output from a LAMMPS log file and plot key diagnostics
(potential energy, temperature, pressure) vs timestep.

Usage:
    python lammpslog_to_csv.py /path/to/*

Expects a thermo_style like:
    thermo_style custom step pe press ke temp lx ly lz pxx pyy pzz spcpu density
"""
import os
import sys
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_thermo_blocks(logfile_dir, columns):
    """
    Extract every thermo table in a LAMMPS log as a separate block
    (e.g. minimization, equilibration, production each get their own).
    Returns a list of lists-of-rows (each row is a list of floats).
    """
    blocks = []
    current_rows = []
    in_table = False

    logfile = os.path.join(logfile_dir, "lammps_run.log")

    with open(logfile) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("Step"):
                in_table = True
                current_rows = []
                continue
            if in_table:
                if stripped.startswith("Loop time"):
                    in_table = False
                    if current_rows:
                        blocks.append(current_rows)
                    continue
                parts = stripped.split()
                if len(parts) == len(columns):
                    try:
                        current_rows.append([float(x) for x in parts])
                    except ValueError:
                        pass  # header re-print or stray non-numeric line
    return blocks


def save_csv(rows, columns, out_csv):
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(rows)
    print(f"saved {out_csv}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python lammpslog_to_csv.py /path/to/*")
        sys.exit(1)

    logfile = sys.argv[1]
    logfile_split = logfile.split("/")
    logfile_dir_name = logfile_split[-1]

    script_dir = os.path.dirname(os.path.abspath(__file__))

    output_location = os.path.join(script_dir, "equil_results")

    save_location_name = os.path.join(output_location, logfile_dir_name)

    # Update this to match your thermo_style custom column order exactly.
    columns = ["Step", "PotEng", "Press", "KinEng", "Temp",
               "Lx", "Ly", "Lz", "Pxx", "Pyy", "Pzz", "SPCPU", "Density"]

    blocks = parse_thermo_blocks(logfile, columns)
    print(f"found {len(blocks)} thermo block(s)")

    for i, rows in enumerate(blocks):
        print(f"  block {i}: {len(rows)} rows, "
              f"steps {rows[0][0]:.0f} to {rows[-1][0]:.0f}")
        save_csv(rows, columns, save_location_name+".csv")