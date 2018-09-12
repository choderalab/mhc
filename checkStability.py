from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

def check_stability(pdbFn, xmlFn):
    pdb = PDBFile(pdbFn)
    with open(xmlFn, 'r') as f:
       system = XmlSerializer.deserialize(f.read())
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()
    simulation.reporters.append(PDBReporter('Complexoutput.pdb', 1000))
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
            potentialEnergy=True, temperature=True))
    simulation.step(10000)

if __name__ == '__main__':
    from sys import argv
    from subroutines import get_opts

    opts = get_opts(argv)
    pdbFn = opts['-pdbfn']
    xmlFn = opts['-xmlfn']
    check_stability(pdbFn, xmlFn)