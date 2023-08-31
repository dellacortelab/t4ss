from openmm.unit import *
from openmm import *
from openmm.app import *
import martini_openmm as martini
from mdtraj.reporters import XTCReporter # This doesn't seem to exist on the internet, but I don't get errors
# from checkpointreporter import CheckpointReporter ## I think it can already access this
from sys import stdout


# This script originates as a copy of run.py from Brenden's OpenMM work on the complex without a membrane.
# My goal is both to compare this run to an equivalent gromacs simulation and prepare to do an alchemical simulation with the whole complex.

def run(epsilon_r):

	platform = Platform.getPlatformByName("CUDA")
	properties = {'Precision': 'mixed', 'UseCpuPme': 'true'}
	# I'm given to believe these properties may speed up the process, which will be critical when running an alchemical simulation.

	conf = GromacsGroFile("../system.gro")
	box_vectors = conf.getPeriodicBoxVectors()

	# get any defines
	# This is vestigial from the copy. I'm not sure if we need defines.
	defines = {}
	try:
		with open("../../defines.txt") as def_file:
			for line in def_file:
				line = line.strip()
				defines[line] = True
	except FileNotFoundError:
		pass

	top = martini.MartiniTopFile(
		"../system.top",
		periodicBoxVectors=box_vectors,
		defines=defines,
		epsilon_r=epsilon_r,
	)
	system = top.create_system(nonbonded_cutoff=1.1 * nanometer)

	integrator = LangevinIntegrator(310 * kelvin,
									10.0 / picosecond,
									20 * femtosecond)
	integrator.setRandomNumberSeed(0)

	simulation = Simulation(top.topology, system, integrator,
							platform, properties)

	simulation.context.setPositions(conf.getPositions())

################################################################################
		### Minimization ###

	simulation.reporters.append(PDBReporter('mini.pdb', 1000))
	simulation.reporters.append(StateDataReporter(stdout, 5000,
													step=True,
													potentialEnergy=True,
													temperature=True,
													volume=True)
								)
	print("Minimizing energy...")
	simulation.minimizeEnergy(maxIterations=5000,tolerance=1.0)

	energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
	print("System minimized at", energies, "\n")

################################################################################
		### NVT equilibration ###

	simulation.context.setVelocitiesToTemperature(310 * kelvin)
	print('Running NVT equilibration...')
	simulation.step(50000) #1ns

################################################################################
		### NPT equilibration ###
	
	system.addForce(MonteCarloBarostat(1 * bar, 310 * kelvin))
	# to update the simulation object to take in account the new system
	simulation.context.reinitialize(True)
	print('Running NPT equilibration...')
	simulation.step(50000) #1ns

	# save the equilibration results to file
	simulation.saveState('equi.state')
	simulation.saveCheckpoint('equi.chk')

################################################################################
		### Production run ###

	# Set up the reporters to report energies every 1000 steps.
	simulation.reporters.append(StateDataReporter("prod.log", 1000,
													step=True,
													potentialEnergy=True,
													totalEnergy=True,
													density=True,
													temperature=True,
													volume=True)
								)
	# save the trajectory in XTC format
	xtc_reporter = XTCReporter('prod.xtc', 10000)
	simulation.reporters.append(xtc_reporter)

	# save the checkpoint every hour real time
	# Steps in first sim = 23M >> 23M / 72 hr ~ 300K / hr (about 3 ns / hr)
	chk_reporter = CheckpointReporter('prod.chk', 300000)
	simulation.reporters.append(chk_reporter)

	# run simulation
	print("Running simulation...")
	simulation.step(50000000) # 1000 ns


run(15)
# why 15? Not sure.
