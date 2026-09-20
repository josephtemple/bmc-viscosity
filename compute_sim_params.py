"""
Recreate Fig. 6d (Joseph et al., Mpipi paper) from LAMMPS fix ave/chunk output.

Expects arg.density.profile, lys.density.profile, rna.density.profile
in the current directory (same format as your fix ave/chunk output).
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


def average_equilibrated(blocks, frac_discard=0.25):
    """Average density over the last (1-frac_discard) fraction of blocks."""
    n = len(blocks)
    keep = blocks[int(n * frac_discard):]
    coords = keep[0][1]
    dens = np.mean([b[2] for b in keep], axis=0)
    return coords, dens


def find_new_box_size(rho_arg, rho_lys, rho_rna, n_arg=50, m_arg = 156.2, m_lys = 128.2, m_rna = 306.2,
                      arg_per_chain = 50, lys_per_chain = 50, rna_per_chain = 10):
    """
    Determine the necessary box size of isotropic simulation (for viscosity) based on dense phase
    densities of arg, lys, rna. Assuming the number of arginine chains is 128 (same as before)

    rho_i = n_i * mass_i / L^3, and L^3 is constant so we find the necessary numbers
    of each species, then solve for L.

    rho given in g/cm^3
    mass given in g/mol, and is the mass of one residue (amino acid/nucleotide). so then multiply
        by the number in a chain to get the mass of the chain
    n is the number of chains
    """
    m_arg_chain = m_arg * arg_per_chain
    m_lys_chain = m_lys * lys_per_chain
    m_rna_chain = m_rna * rna_per_chain

    n_lys = int(n_arg * (m_arg_chain / rho_arg) * (rho_lys / (m_lys_chain)) )
    n_rna = int(n_arg * (m_arg_chain / rho_arg) * (rho_rna / (m_rna_chain)) )

    L = (n_arg * m_arg_chain * CONV / rho_arg)**(1/3)

    return L, {"arg": n_arg, "lys": n_lys, "rna": n_rna}



def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    profile_dir = os.path.join(script_dir, "josephetal22")
    species_files = {
        "rna": os.path.join(profile_dir, "rna.density.profile"),
        "arg": os.path.join(profile_dir, "arg.density.profile"),
        "lys": os.path.join(profile_dir, "lys.density.profile"),
    }

    results = {}
    for name, fn in species_files.items():
        blocks = parse_chunk_file(fn)
        print(f"{name}: {len(blocks)} blocks parsed from {fn} "
              f"(first step {blocks[0][0]}, last step {blocks[-1][0]})")
        coords, dens = average_equilibrated(blocks, frac_discard=0.5)
        results[name] = (coords, dens * CONV)

    # center around peak arginine (as per paper)
    arg_coords, arg_dens = results["arg"]
    binsize = arg_coords[1] - arg_coords[0]
    Lz = arg_coords[-1] + binsize  # approx box length along z
    peak_coord = arg_coords[np.argmax(arg_dens)]
    shift = Lz / 2 - peak_coord

    core_density = {}
    for name in ["rna", "arg", "lys"]:
        coords, dens = results[name]
        x = (coords + shift) % Lz
        order = np.argsort(x)

        # extract densities in core
        pct_range = 0.05 # half width of flat core, as a percentage of total box length. so 0.05 means dense phase is 10% of box
        core_mask = (x[order] > (0.5 - pct_range)*Lz) & (x[order] < (0.5 + pct_range)*Lz)
        ordered_dens = dens[order]
        core_dens = np.mean(ordered_dens[core_mask])

        core_density[name] = core_dens
        print(name, core_density[name])

    new_box_length, counts = find_new_box_size(core_density["arg"], core_density["lys"], core_density["rna"])
    
    print(counts["arg"], "arginine chains")
    print(counts["lys"], "lysine chains")
    print(counts["rna"], "rna chains")
    print("L =", new_box_length)

    # write to file
    visc_path = os.path.dirname(script_dir)
    write_to = os.path.join(visc_path, "sim_params.txt")
    with open(write_to, "w") as file:
        file.write(f"arg density (g/cm^3): {core_density['arg']} \n")
        file.write(f"lys density (g/cm^3): {core_density['lys']} \n")
        file.write(f"rna density (g/cm^3): {core_density['rna']} \n")
        file.write(f"sim box length (ang): {new_box_length}\n")
        file.write(f"num arg needed:       {counts['arg']}\n")
        file.write(f"num lys needed:       {counts['lys']}\n")
        file.write(f"num rna needed:       {counts['rna']}\n")



if __name__ == "__main__":
    main()