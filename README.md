# hladr1-tpi
Goal: set up complex and solvent files {xml, pdb} for simulation in openmm
1. From RCSB PDB file, solvate, minimize energy, and write out serialized xml and pdb files.
    - as much resolved density as possible
    - check for missing loops
    - check comments
    - for both peptide alone and peptide with MHC
2. Basic simulations in openmm to ensure stability
