
from simtk import openmm, unit
from simtk.openmm.app import PDBFile
from openff.toolkit.topology import Molecule
from simtk.openmm.app import AmberPrmtopFile

import sys
import time

pdbfile = PDBFile(sys.argv[1])
prmtop = AmberPrmtopFile(sys.argv[2])

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.2*nanometer,constraints=HBonds)

# Propagate the System with Langevin dynamics.
time_step = 2*unit.femtoseconds  # simulation timestep
temperature = 300*unit.kelvin  # simulation temperature
friction = 1/unit.picosecond  # collision rate
integrator = openmm.LangevinIntegrator(temperature, friction, time_step)

# Length of the simulation.
num_steps = 10000  # number of integration steps to run

# Logging options.
trj_freq = 1  # number of steps per written trajectory frame
data_freq = 1  # number of steps per written simulation statistics

# Set up an OpenMM simulation.
simulation = openmm.app.Simulation(prmtop.topology, system, integrator)

# Set the initial positions.
positions = pdbfile.getPositions() 
simulation.context.setPositions(positions)

# Randomize the velocities from a Boltzmann distribution at a given temperature.
simulation.context.setVelocitiesToTemperature(temperature)

# Configure the information in the output files.
pdb_reporter = openmm.app.PDBReporter('trajectory.pdb', trj_freq)
state_data_reporter = openmm.app.StateDataReporter('data.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True)
simulation.reporters.append(pdb_reporter)
simulation.reporters.append(state_data_reporter)

print("Starting simulation")
start = time.process_time()

# Run the simulation
simulation.step(num_steps)

end = time.process_time()
print("Elapsed time %.2f seconds" % (end-start))
print("Done!")

