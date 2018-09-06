#!/usr/bin/env python

"""
Retrieve a PDB file from the RCSB, solvate it, and minimize.

"""

################################################################################
# IMPORTS
################################################################################

from simtk import unit
from simtk import openmm
from simtk.openmm import app
from pdbfixer import PDBFixer

################################################################################
# OPTIONS
################################################################################

pdbid = '1klu'  # PDB ID to retrieve
chain_ids_to_keep = ['A', 'B', 'C']  # chains to keep
if len(chain_ids_to_keep) > 1:
    fnPrefix = pdbid + 'Complex'
else:
    fnPrefix = pdbid + 'Ligand'
pH = 7.0 # pH
forcefields_to_use = ['amber99sbildn.xml', 'tip3p.xml'] # list of forcefields to use in parameterization
padding = 11.0 * unit.angstroms # padding to use for adding solvent
nonbondedMethod = app.PME # nonbonded method
constraints = app.HBonds # bonds to be constrained
keepWater = False # keep crystal water

################################################################################
# SUBROUTINES
################################################################################

def write_file(filename, contents):
    outfile = open(filename, 'w')
    outfile.write(contents)
    outfile.close()

################################################################################
# SET UP SYSTEM
################################################################################

# Load forcefield.
forcefield = app.ForceField(*forcefields_to_use)

# Retrieve structure from PDB.
print('Retrieving %s from PDB...' % pdbid)
fixer = PDBFixer(pdbid=pdbid)

# Build a list of chains to remove.
print('Removing all chains but %s' % chain_ids_to_keep)
all_chains = list(fixer.topology.chains())
chain_id_list = [c.id for c in fixer.topology.chains()]
chain_ids_to_remove = set(chain_id_list) - set(chain_ids_to_keep)
fixer.removeChains(chainIds=chain_ids_to_remove)

# Find missing residues.
print('Finding missing residues...')
fixer.findMissingResidues()

# Replace nonstandard residues.
print('Replacing nonstandard residues...')
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()

# Add missing atoms.
print('Adding missing atoms...')
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Remove heterogens.
print('Removing heterogens...')
fixer.removeHeterogens(keepWater=keepWater)

# Add missing hydrogens.
print('Adding missing hydrogens appropriate for pH %s' % pH)
fixer.addMissingHydrogens(pH)

# Write PDB file.
output_filename = '%s-pdbfixer-nosolvent.pdb' % fnPrefix
print('Writing PDB file to "%s"...' % output_filename)
app.PDBFile.writeFile(fixer.topology, fixer.positions, open(output_filename, 'w'))

if nonbondedMethod in [app.PME, app.CutoffPeriodic, app.Ewald]:
    # Add solvent.
    print('Adding solvent...')
    fixer.addSolvent(padding=padding)

# Write PDB file.
output_filename = '%s-pdbfixer.pdb' % fnPrefix
print('Writing PDB file to "%s"...' % output_filename)
app.PDBFile.writeFile(fixer.topology, fixer.positions, open(output_filename, 'w'))

# Create OpenMM System.
print('Creating OpenMM system...')
system = forcefield.createSystem(fixer.topology, nonbondedMethod=nonbondedMethod,
                                 constraints=constraints, rigidWater=True, removeCMMotion=False)

# Minimimze to update positions.
print('Minimizing...')
integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
context = openmm.Context(system, integrator)
context.setPositions(fixer.positions)
openmm.LocalEnergyMinimizer.minimize(context)
state = context.getState(getPositions=True)
fixer.positions = state.getPositions()

# Write final coordinates.
output_filename = '%s-minimized.pdb' % fnPrefix
print('Writing PDB file to "%s"...' % output_filename)
app.PDBFile.writeFile(fixer.topology, fixer.positions, open(output_filename, 'w'))

# Serialize final coordinates.
print('Serializing to XML...')
system_filename = fnPrefix + 'System.xml'
integrator_filename = fnPrefix + 'Integrator.xml'
state_filename = fnPrefix + 'State.xml'
write_file(system_filename, openmm.XmlSerializer.serialize(system))
write_file(integrator_filename, openmm.XmlSerializer.serialize(integrator))
state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True,
                         getParameters=True, enforcePeriodicBox=True)
write_file(state_filename, openmm.XmlSerializer.serialize(state))

