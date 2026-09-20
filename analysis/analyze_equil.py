"""
Parse thermo output from a LAMMPS log file and plot key diagnostics
(potential energy, temperature, pressure) vs timestep.

Usage:
    python analyze_equil.py /path/to/lammps_run.log 

Expects a thermo_style like:
    thermo_style custom step pe press ke temp lx ly lz pxx pyy pzz spcpu density
"""
import sys
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_thermo_blocks(logfile, columns):
    """
    Extract every thermo table in a LAMMPS log as a separate block
    (e.g. minimization, equilibration, production each get their own).
    Returns a list of lists-of-rows (each row is a list of floats).
    """
    blocks = []
    current_rows = []
    in_table = False

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


def plot_block(rows, columns, out_png):
    idx = {name: i for i, name in enumerate(columns)}
    steps = [r[idx["Step"]] for r in rows]

    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)

    axes[0].plot(steps, [r[idx["PotEng"]] for r in rows], color="tab:blue")
    axes[0].set_ylabel("Potential Energy")
    axes[0].set_title("Thermo diagnostics")

    axes[1].plot(steps, [r[idx["Temp"]] for r in rows], color="tab:red")
    axes[1].set_ylabel("Temperature (K)")

    axes[2].plot(steps, [r[idx["Press"]] for r in rows], color="tab:green", label="Press (avg)")
    if "Pxx" in idx:
        axes[2].plot(steps, [r[idx["Pxx"]] for r in rows], alpha=0.4, label="Pxx")
        axes[2].plot(steps, [r[idx["Pyy"]] for r in rows], alpha=0.4, label="Pyy")
        axes[2].plot(steps, [r[idx["Pzz"]] for r in rows], alpha=0.4, label="Pzz")
    axes[2].set_ylabel("Pressure (atm)")
    axes[2].set_xlabel("Timestep")
    axes[2].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    print(f"saved {out_png}")


def save_csv(rows, columns, out_csv):
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(rows)
    print(f"saved {out_csv}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python run_equil.py /path/to/lammps_run.log")
        sys.exit(1)

    logfile = sys.argv[1]

    # Update this to match your thermo_style custom column order exactly.
    columns = ["Step", "PotEng", "Press", "KinEng", "Temp",
               "Lx", "Ly", "Lz", "Pxx", "Pyy", "Pzz", "SPCPU", "Density"]

    blocks = parse_thermo_blocks(logfile, columns)
    print(f"found {len(blocks)} thermo block(s)")

    for i, rows in enumerate(blocks):
        print(f"  block {i}: {len(rows)} rows, "
              f"steps {rows[0][0]:.0f} to {rows[-1][0]:.0f}")
        save_csv(rows, columns, f"thermo_block{i+1}.csv")
        plot_block(rows, columns, f"thermo_block{i+1}.png")