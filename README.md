# BMC-Viscosity
LAMMPS simulations of dense-phase biomolecular condensates measuring material properties
Python scripts to generate input files for dense-phase simulations and analysis as we go.

Research currently in progress in the GrandPre lab at Washington University in St. Louis.

As is, we borrow results of the direct-coexistence simulation from "Physics-driven coarse-grained model for biomolecular phase separation with near-quantitative accuracy" by Jerelle Joseph et al 2022 to measure the densities of a PolyR-PolyK-PolyU system in the dense phase. From there, my work includes Python scripts to make the rest of the simulation pipeline (isotropic dense-phase simulation construction, equilibration, into MD simulation measured mean square displacements for diffusion constant and stress tensor for viscosity) as seamless as it can be.
