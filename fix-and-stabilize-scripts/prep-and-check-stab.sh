#!/bin/bash

# set variables
PDB_ID='1klu'
COMPLEX_CHAINS='A,B,C'
LIGAND_CHAINS='C'

# prepare PDB files
echo 'preparing complex PDB files'
python prepPdbFile.py -pdbid $PDB_ID -chains $COMPLEX_CHAINS -id complex
echo 'preparing ligand PDB files'
python prepPdbFile.py -pdbid $PDB_ID -chains $LIGAND_CHAINS -id ligand

# check stability
LIGAND_PDB_FN=$PDB_ID'-ligand-minimized.pdb'
LIGAND_XML_FN=$PDB_ID'-ligand-system.xml'
COMPLEX_PDB_FN=$PDB_ID'-complex-minimized.pdb'
COMPLEX_XML_FN=$PDB_ID'-complex-system.xml'

echo 'checking stability of ligand'
python checkStability.py -pdbfn $LIGAND_PDB_FN -xmlfn $LIGAND_XML_FN > ligand_stability.out
echo 'checking stability of complex'
python checkStability.py -pdbfn $COMPLEX_PDB_FN -xmlfn $COMPLEX_XML_FN > complex_stability.out
