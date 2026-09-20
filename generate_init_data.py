"""
Generate initial data for isotropic simulation box using Mpipi parameters, saving to
./dense_phase.dat


# ATOM LINE is: atom-id mol-id atom-type charge x-coord y-coord z-coord ix iy iz

where atom-id a unique identifier to that atom, mol-id is a unique identifier to that chain
atom type is the number (5,3,44)->(arg,ls,rna), charge is the charge, coords are 
the coords as expected, and ix iy iz are "image flags" that tell you how many times the coordinate
had to be wrapped aroud the given box of size L.

i_n = floor(coord / L)   # integer, per-axis
wrapped coord = coord - i_n * L


# BOND LINE is: bond-id bond-type atom1 atom2
bond-id is a unique identifier to that bond, bond-type is either {1,2} where 1 is used for
polypeptide chains and 2 is used for rna chains. atom1 and atom2 are the atom-id of the two
molecules with the bond.

"""
import numpy as np
import os

def random_unit_vector():
    r = np.array(np.random.normal(size=3))
    r /= np.linalg.norm(r)
    return r
    
# construct atom coordinates. if its the first one in the polymer, place randomly in the box
# place subsequent atoms by stepping by the bond length in a random direction
def atom_coords(L, bond_length, prev_pos = None):
    # if bead is the first in its chain, choose random place in box
    if prev_pos is None:
        pos = np.random.uniform(0, L, size=3)
    else:
        pos = prev_pos + bond_length * random_unit_vector()

    image = np.floor(pos / L).astype(int)
    wrapped_pos = pos - image * L        

    return wrapped_pos, image

# eventually replace this to take a commandline argument of name for dense_phase.dat
def main():
    # replace at some point with file read (split each line at spaces and grab the [-1] element)
    L = 106.89119293174846

    num_arg = 50
    num_lys = 0
    num_rna = 280

    arg_per_chain = 50
    lys_per_chain = 50
    rna_per_chain = 10

    # dictionaries of necessary things
    atom_type_dict = {"arg": 5, "lys": 3, "rna": 44} # in mpipi, arginine is labeled 5, lysine is 3, and uracil is 44

    counts_dict = {"arg": num_arg, "lys": num_lys, "rna": num_rna}
    per_chain_counts_dict = {"arg": arg_per_chain, "lys": lys_per_chain, "rna": rna_per_chain}

    total_num_atoms = 0
    total_num_bonds = 0
    for species in ['arg','lys','rna']:
        total_num_atoms += counts_dict[species] * per_chain_counts_dict[species]
        total_num_bonds += counts_dict[species] * (per_chain_counts_dict[species] - 1)   # chain with n beads has n-1 bonds

    bond_type_dict = {"arg": 1, "lys": 1, "rna": 2}
    bond_lengths_dict = {"arg": 3.81, "lys": 3.81, "rna": 5.0} # as defined in the harmonic potential of mpipi
    charge_dict = {"arg": 0.75, "lys": 0.75, "rna": -0.75}


    # write to .dat file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(script_dir, "dense_phase.dat")
    with open(input_file, "w") as file:
        file.write(f"LAMMPS data file for condensate dense phase simulation\n\n")
        file.write(f"{total_num_atoms} atoms\n")
        file.write("44 atom types\n")
        file.write(f"{total_num_bonds} bonds\n")
        file.write("2 bond types\n\n")
        for n in ['x','y','z']:
            file.write(f"0.0 {L} {n}lo {n}hi\n")
        file.write("\nMasses\n\n")

        masses_file = os.path.join(script_dir, "mpipi_masses.dat")
        with open(masses_file, "r") as mass_file:
            for line in mass_file:
                file.write(f"{line.strip()}\n") 

        file.write("\nAtoms\n\n")

        atom_id = 1
        mol_id = 1

        # loop thru different species
        # in type, loop thru chains
        # in chain, loop thru atoms
        for species in ['arg','lys','rna']:
            atom_type = atom_type_dict[species]
            charge = charge_dict[species]
            for chain in range(counts_dict[species]):
                # do stuff
                first_in_chain = True
                for atom in range(per_chain_counts_dict[species]):
                    if first_in_chain:
                        wrapped_pos, img = atom_coords(L=L, bond_length=bond_lengths_dict[species])
                        first_in_chain = False
                    else:
                        # images need true position to be computed consistently
                        wrapped_pos, img = atom_coords(L=L, bond_length=bond_lengths_dict[species], prev_pos=true_pos)

                    # positions wrapped within box coordinates
                    x, y, z = wrapped_pos
                    ix, iy, iz = img
                    file.write(f"{atom_id} {mol_id} {atom_type} {charge} {x} {y} {z} {ix} {iy} {iz}\n")

                    # to be used in next loop
                    true_pos = wrapped_pos + img * L   # reconstruct the real, unwrapped position
                    
                    atom_id += 1
                mol_id += 1

        # loop thru
        file.write("\nBonds\n\n")
        
        # loop thru different species
        # in type, loop thru chains
        # in chain, loop thru atoms
        atom_id = 1
        bond_id = 1

        for species in ['arg','lys','rna']:
            bond_type = bond_type_dict[species]
            for chain in range(counts_dict[species]):
                chain_start = atom_id
                for i in range(per_chain_counts_dict[species] - 1):
                    atom1 = chain_start + i
                    atom2 = chain_start + i + 1
                    file.write(f"{bond_id} {bond_type} {atom1} {atom2}\n")
                    bond_id += 1
                atom_id += per_chain_counts_dict[species]
        


if __name__ == "__main__":
    main()