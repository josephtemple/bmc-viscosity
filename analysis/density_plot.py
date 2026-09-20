"""
Recreate Fig. 6d (Joseph et al., Mpipi paper) from LAMMPS fix ave/chunk output.

Expects arg.density.profile, lys.density.profile, rna.density.profile
in the ../josephetal22 directory.
"""
import numpy as np
import matplotlib.pyplot as plt
import os

# g/mol per Angstrom^3  ->  g/cm^3
CONV = 1e24 / 6.02214076e23

def parse_chunk_file(path):
    """Parse a LAMMPS fix ave/chunk output file.

    Returns a list of (timestep, coords, density_mass) tuples, one per
    output block. coords and density_mass are numpy arrays of length
    nchunks.
    """
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f if l.strip() and not l.strip().startswith("#")]

    blocks = []
    i = 0
    while i < len(lines):
        header = lines[i].split()
        step, nchunks, _total = int(header[0]), int(header[1]), int(header[2])
        coords = np.zeros(nchunks)
        density_mass = np.zeros(nchunks)
        for j in range(nchunks):
            parts = lines[i + 1 + j].split()
            # columns: chunk_id coord1 ncount density/mass density/number
            coords[j] = float(parts[1])
            density_mass[j] = float(parts[3])
        blocks.append((step, coords, density_mass))
        i += nchunks + 1
    return blocks


def average_equilibrated(blocks, frac_discard=0.0):
    """Average density over the last (1-frac_discard) fraction of blocks."""
    n = len(blocks)
    keep = blocks[int(n * frac_discard):]
    coords = keep[0][1]
    dens = np.mean([b[2] for b in keep], axis=0)
    return coords, dens


def main():
    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "stix"
    })
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)
    profiles_dir = os.path.join(parent_dir, "josephetal22")
    species_files = {
        "rna": os.path.join(profiles_dir, "rna.density.profile"),
        "arg": os.path.join(profiles_dir, "arg.density.profile"),
        "lys": os.path.join(profiles_dir, "lys.density.profile"),
    }
    colors = {"rna": "gold", "arg": "magenta", "lys": "green"}

    results = {}
    for name, fn in species_files.items():
        blocks = parse_chunk_file(fn)
        print(f"{name}: {len(blocks)} blocks parsed from {fn} "
              f"(first step {blocks[0][0]}, last step {blocks[-1][0]})")
        coords, dens = average_equilibrated(blocks, frac_discard=0.5)
        results[name] = (coords, dens * CONV)

    # Recenter on the arg-rich core (per the paper, arg sits at the droplet
    # centre) so the condensate lines up in the middle of the plot.
    arg_coords, arg_dens = results["arg"]
    binsize = arg_coords[1] - arg_coords[0]
    Lz = arg_coords[-1] + binsize  # approx box length along z
    peak_coord = arg_coords[np.argmax(arg_dens)]
    shift = Lz / 2 - peak_coord

    img_file = os.path.join(script_dir, "density_plot.png")
    plt.figure(figsize=(6, 4))
    for name in ["rna", "arg", "lys"]:
        coords, dens = results[name]
        x = (coords + shift) % Lz
        order = np.argsort(x)
        plt.plot(x[order] / Lz, dens[order], label=name, color=colors[name])

    plt.xlabel("Box long axis, L (normalized)")
    plt.ylabel(r"Density / g cm$^{-3}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(img_file, dpi=300)
    print("Saved density_plot.png")


if __name__ == "__main__":
    main()