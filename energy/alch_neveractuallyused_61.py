from openmm.unit import *
from openmm import *
from openmm.app import *
import martini_openmm as martini
from mdtraj.reporters import XTCReporter
from sys import stdout
from sys import argv
import copy
import numpy as np

# This originates as a copy of alch_3.py. The naming refers to the number of lambdas.
# The goal is to double the distribution's sampling. I wanted to try complementary intercalation of the lambda values, but that wouldn't
# be compatible with the energy matrix that gets assembled further down.

### Establish the lambdas
cindex = 0
vindex = 1
# Length is currently 61
lambda_twople_list = np.array([
        (1,1),
        (0.95,1.00),
        (0.9,1.0),
        (0.85,1.00),
        (0.8,1.0),
        (0.75,1.00),
        (0.7,1.0),
        (0.65,1.00),
        (0.6,1.0),
        (0.55,1.00),
        (0.5,1.0),
        (0.45,1.00),
        (0.4,1.0),
        (0.35,1.00),
        (0.3,1.0),
        (0.25,1.00),
        (0.2,1.0),
        (0.15,1.00),
        (0.1,1.0),
        (0.05,1.00),
        (0,1),
        (0.000,0.975),
        (0,0.95),
        (0.000,0.925),
        (0,0.90),
        (0.000,0.875),
        (0,0.85),
        (0.000,0.825),
        (0,0.80),
        (0.000,0.775),
        (0,0.75),
        (0.000,0.725),
        (0,0.70),
        (0.000,0.675),
        (0,0.65),
        (0.000,0.625),
        (0,0.60),
        (0.000,0.575),
        (0,0.55),
        (0.000,0.525),
        (0,0.50),
        (0.000,0.475),
        (0,0.45),
        (0.000,0.425),
        (0,0.40),
        (0.000,0.375),
        (0,0.35),
        (0.000,0.325),
        (0,0.30),
        (0.000,0.275),
        (0,0.25),
        (0.000,0.225),
        (0,0.20),
        (0.000,0.175),
        (0,0.15),
        (0.000,0.125),
        (0,0.10),
        (0.000,0.075),
        (0,0.05),
        (0.000,0.025),
        (0,0)
])

def run(epsilon_r):


        ### Tuple option
        if len(argv) < 2:
                print("Don't forget the lambda argument")
                return
        place = int(argv[1])
        CLAMBDA = lambda_twople_list[place][cindex]
        VLAMBDA = lambda_twople_list[place][vindex]
        ###

        print("Formatting")

        platform = Platform.getPlatformByName("CUDA")
        properties = {'Precision': 'mixed'}
        # Double precision is slower than mixed. Mixed doesn't appear to cause problems for mbar.
        # For more info on properties: http://docs.openmm.org/latest/userguide/library/04_platform_specifics.html?highlight=mixed

        conf = GromacsGroFile("../system.gro")
        box_vectors = conf.getPeriodicBoxVectors()

        # defines deleted

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

        simulation = Simulation(top.topology, system, integrator, platform, properties)
        print("Force calculations")
        # Retrieve the NonbondedForce
        forces = { force.__class__.__name__ : force for force in system.getForces() }
        nbforce = forces['CustomNonbondedForce']
        print("Deep copying...")
        innerforce = copy.deepcopy(nbforce)


        # Add a CustomNonbondedForce to handle only alchemically-modified interactions

        # Make two sets of particles, one that contains just the particle we will alchemically annihilate
        # and the other which contains all the other particles.
        # Question: the range makes a set including 0-4974. Where are those used? Consulting the GromacsGroFile? After verifying, it seems that OpenMM treats each
        # particle row as the index.
        alchemical_particles = set(range(4975))
        chemical_particles = set(range(system.getNumParticles())) - alchemical_particles

        nonbonded_cutoff=1.1 * unit.nanometer
        # Define the energy function for the CustomNonbondedForce
        # when lambda is 1.0 it is a normal LJ potential, when lambda is 0.0 the interaction vanishes
        # ADAPTED FROM MARTINI.PY -- MARTINI_OPENMM
        # Custom force to handle martini nonbonded terms
        # nbfix-like terms
        # clambda will be the coulombic lambda
        # vlambda will be the vdw lambda
        ## vlambda can hit ES, since coulombic interactions will always be turned off first
        alch_force = CustomNonbondedForce(
                "step(rcut-r)*(vlambda*(LJ-corr) + clambda*ES);"
                "LJ = C12(type1, type2) / reff_sterics^12 - C6(type1, type2) / reff_sterics^6;"
                "reff_sterics = sigma*(0.5*(1.0-vlambda) + (r/sigma)^6)^(1/6);"
                "corr = C12(type1, type2) / reff_cut^12 - C6(type1, type2) / reff_cut^6;"
                "reff_cut = sigma*(0.5*(1.0-vlambda) + (rcut/sigma)^6)^(1/6);"
                "sigma = (C12(type1, type2)/C6(type1, type2))^(1/6);"
                "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
                "crf = 1 / rcut + krf * rcut^2;"
                "krf = 1 / (2 * rcut^3);"
                f"epsilon_r = {epsilon_r};"
                "f = 138.935458;"
                f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
                )
        alch_force.addPerParticleParameter("type")
        alch_force.addPerParticleParameter("q")
        alch_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        alch_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))

        # Add lambda as a parameter we can change during the simulation
        # This is before minimization.
        alch_force.addGlobalParameter('clambda', 1.0)
        alch_force.addGlobalParameter('vlambda', 1.0)

        # set the values of sigma and epsilon by copying them from the existing NonBondedForce
        for index in range(system.getNumParticles()):
                [type, charge] = nbforce.getParticleParameters(index)
                alch_force.addParticle([type, charge])

        # Set the custom force to occur between just the alchemical particle and the other particles
        alch_force.addInteractionGroup(alchemical_particles, chemical_particles)
        system.addForce(alch_force)
        innerforce.addInteractionGroup(alchemical_particles, alchemical_particles)
        system.addForce(innerforce)
        nbforce.addInteractionGroup(chemical_particles, chemical_particles)

        print("Alchemical force added")
################################################################################
        # THE FOLLOWING IS MAINLY DRAWN FROM MARTINI.PY -- MARTINI_OPENMM
        # build a lookup table mapping atom types into integer indices
        # that will later be used to lookup LJ combinations
        used_atom_types = set()
        for molecule_name, _ in top._molecules:
                molecule_type = top._moleculeTypes[molecule_name]
                for atom in molecule_type.atoms:
                        used_atom_types.add(atom[1])
        atom_type_map = {k: i for i, k in enumerate(sorted(used_atom_types))}

        # now we need to setup our tables of C6 and C12
        n_types = len(atom_type_map)
        C6 = []
        C12 = []
        for type_i in atom_type_map:
                for type_j in atom_type_map:
                        type_i_sorted, type_j_sorted = sorted([type_i, type_j])
                        if (type_i_sorted, type_j_sorted) in top._nonbond_types:
                                # Parameters were specified for this pair of types,
                                # so use them.
                                params = top._nonbond_types[(type_i_sorted, type_j_sorted)]
                                if top._use_sigma_eps:
                                        sigma = float(params[3])
                                        eps = float(params[4])
                                        c6 = 4 * eps * sigma ** 6
                                        c12 = 4 * eps * sigma ** 12
                                else:
                                        c6 = float(params[3])
                                        c12 = float(params[4])
                        else:
                                # Parameters were not specified for this pair type,
                                # so calculate using combination rules.
                                params_i = top._atom_types[type_i]
                                v_i = float(params_i[6])
                                w_i = float(params_i[7])
                                params_j = top._atom_types[type_j]
                                v_j = float(params_j[6])
                                w_j = float(params_j[7])

                                if top._use_sigma_eps:
                                        sigma = 0.5 * (v_i + v_j)
                                        eps = math.sqrt(w_i * w_j)
                                        c6 = 4 * eps * sigma ** 6
                                        c12 = 4 * eps * sigma ** 12
                                else:
                                        c6 = math.sqrt(v_i * v_j)
                                        c12 = math.sqrt(w_i * w_j)

                        C6.append(c6)
                        C12.append(c12)
        alch_force.addTabulatedFunction(
                "C6", Discrete2DFunction(n_types, n_types, C6)
        )
        alch_force.addTabulatedFunction(
                "C12", Discrete2DFunction(n_types, n_types, C12)
        )
        print("Alchemical force complete with 2d functions")

        
        print("Simulation initialized")
        simulation.context.setPositions(conf.getPositions())

################################################################################
                ### Minimization ###

        # simulation.reporters.append(PDBReporter('mini.pdb', 1000)) # This seems to record the production run.
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
        
        simulation.context.reinitialize(True)
        simulation.context.setParameter('clambda', CLAMBDA)
        simulation.context.setParameter('vlambda', VLAMBDA)

        simulation.context.setVelocitiesToTemperature(310 * kelvin)
        print('Running NVT equilibration...')
        simulation.step(5000) #100ps

################################################################################
                ### NPT equilibration ###

        system.addForce(MonteCarloBarostat(1 * bar, 310 * kelvin))
        # to update the simulation object to take in account the new system
        simulation.context.reinitialize(True)
        print('Running NPT equilibration...')
        simulation.step(5000) #100ps

        # save the equilibration results to file
        simulation.saveState('equi.state')
        simulation.saveCheckpoint('equi.chk')

################################################################################
                ### Production run ###

        # Collect data
        xtc_reporter = XTCReporter('prod.xtc', 10000)
        simulation.reporters.append(xtc_reporter)

        ###
        # Thoughts for the following section.
        # Each run of alch.py will be for a particular lambda index.
        # So each run of alch.py will be responsible for one k. 
        # Do we still do multiple iterations?
        # May as well assume so until further notice.
        ###

        nsteps = 2500
        niterations = 100
        nstates = len(lambda_twople_list)

        # u_kln
        # Each lambda gets its own u_kln. Needs to be concatenated with the rest before passing it to pymbar.
        u_kln = np.zeros([nstates,nstates,niterations], np.float64)

        kT = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB * integrator.getTemperature()

        # for k in range(nstates):
        k = place
        for iteration in range(niterations):
                print('state %5d iteration %5d / %5d' % (k, iteration, niterations))
        
                # Run some dynamics
                simulation.step(nsteps)
                # Compute energies at all alchemical states
                for l in range(nstates):
                        simulation.context.setParameter('clambda', lambda_twople_list[l][cindex])
                        simulation.context.setParameter('vlambda', lambda_twople_list[l][vindex])
                        u_kln[k,l,iteration] = simulation.context.getState(getEnergy=True).getPotentialEnergy() / kT
                
                # Reset alchemical state
                simulation.context.setParameter('clambda', CLAMBDA)
                simulation.context.setParameter('vlambda', VLAMBDA)

        newarr = u_kln.reshape(u_kln.shape[0], (u_kln.shape[1]*u_kln.shape[2]))
        np.savetxt("u_kln.txt",newarr)

run(15)
