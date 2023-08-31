from openmm.unit import *
from openmm import *
from openmm.app import *
import martini_openmm as martini
from mdtraj.reporters import XTCReporter
# from checkpointreporter import CheckpointReporter ## I think it can already access this
from sys import stdout

THIS_ROUND = 3

# This script originates as a copy of run.py from Brenden's OpenMM work on the complex without a membrane.
# My goal is both to compare this run to an equivalent gromacs simulation and prepare to do an alchemical simulation with the whole complex.
# I think.

def run(epsilon_r):

	platform = Platform.getPlatformByName("CUDA")
	properties = {'Precision': 'mixed', 'UseCpuPme': 'true'}
	# I'm given to believe these properties may speed up the process, which will be critical when running an alchemical simulation.

	conf = GromacsGroFile("../system.gro")
	box_vectors = conf.getPeriodicBoxVectors()

	# no defines

	top = martini.MartiniTopFile(
		"../system.top",
		periodicBoxVectors=box_vectors,
		epsilon_r=epsilon_r,
	)
	system = top.create_system(nonbonded_cutoff=1.1 * nanometer)

	integrator = LangevinIntegrator(310 * kelvin,
									10.0 / picosecond,
									20 * femtosecond)
	integrator.setRandomNumberSeed(0)

	simulation = Simulation(top.topology, system, integrator,
							platform, properties)

	# simulation.context.setPositions(conf.getPositions())

################################################################################

################################################################################
		### Production run ###

	# Set up the reporters to report energies every 1000 steps.
	simulation.reporters.append(StateDataReporter("prod.log", 1000,
													step=True,
													potentialEnergy=True,
													totalEnergy=True,
													density=True,
													temperature=True,
													volume=True,
													append=True)
								)
	# save the trajectory in XTC format
	xtc_reporter = XTCReporter('prod.xtc', 10000, append=True)
	simulation.reporters.append(xtc_reporter)

	# save the checkpoint every hour real time
	# Steps in first sim = 23M / 72 hr ~ 300K / hr (about 3 ns / hr)
	# chk_reporter = CheckpointReporter('prod.chk', 146000) ### This is how it was for the ones that crashed on storage
	chk_reporter = CheckpointReporter(f'prod{THIS_ROUND}.chk', 300000) ### Just in case there's a create/load order issue
	simulation.reporters.append(chk_reporter)

	# run simulation
	print("Running simulation...")
	simulation.loadCheckpoint(f'prod{THIS_ROUND-1}.chk')
	simulation.step(50000000) # 1000 ns


run(15)
# why 15? Not sure.
