from openmm.unit import *
from openmm import *
from openmm.app import *
import martini_openmm as martini
from mdtraj.reporters import XTCReporter
from sys import stdout
from sys import argv
import copy
import numpy as np

# This originates as a copy of alch.py
# The goal with this version is to adjust the script to be more in line with the Gromacs tutorial for protein-ligand binding.

### OPTION
cindex = 0
vindex = 1
# Length is currently 31
lambda_twople_list = np.array([
        (1,1),
        (0.9,1.0),
        (0.8,1.0),
        (0.7,1.0),
        (0.6,1.0),
        (0.5,1.0),
        (0.4,1.0),
        (0.3,1.0),
        (0.2,1.0),
        (0.1,1.0),
        (0,1),
        (0,0.95),
        (0,0.90),
        (0,0.85),
        (0,0.80),
        (0,0.75),
        (0,0.70),
        (0,0.65),
        (0,0.60),
        (0,0.55),
        (0,0.50),
        (0,0.45),
        (0,0.40),
        (0,0.35),
        (0,0.30),
        (0,0.25),
        (0,0.20),
        (0,0.15),
        (0,0.10),
        (0,0.05),
        (0,0)
])
###

def run(epsilon_r):

        ###
        # This chunk here is one option. See line 232 for a bit of a challenge (the innermost for loop question).
        # Maybe bring in an array of twoples, describing the clambda and vlambda, along with an index to describe which simulation is being executed.
        # Then the parameter setting can be done that way.
        # But why is that loop even there?
        # Well, it's part of how we pass info to pymbar.
        # if len(argv) < 3:
        #         print("Don't forget the arguments")
        #         return
        # CLAMBDA = argv[1]
        # VLAMBDA = argv[2]
        ###

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
        properties = {'Precision': 'double'}
        # Double precision is slower than mixed
        # For more info on properties: http://docs.openmm.org/latest/userguide/library/04_platform_specifics.html?highlight=mixed

        conf = GromacsGroFile("system.gro")
        box_vectors = conf.getPeriodicBoxVectors()

        # get any defines
        defines = {}
        try:
                with open("defines.txt") as def_file:
                        for line in def_file:
                                line = line.strip()
                                defines[line] = True
        except FileNotFoundError:
                pass

        top = martini.MartiniTopFile(
                "system.top",
                periodicBoxVectors=box_vectors,
                defines=defines,
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
        # Question: the range makes a set including 0-4974. Where are those used? Consulting the GromacsGroFile?
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
                "step(rcut-r)*(LJ - corr + ES)*vlambda;"
                "LJ = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
                "corr = (C12(type1, type2) / rcut^12 - C6(type1, type2) / rcut^6);"
                "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf) * clambda;"
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
        # FIXME: This is before minimization. Probably have to go back to "reinitialize" after minimization.
        alch_force.addGlobalParameter('clambda', CLAMBDA)
        alch_force.addGlobalParameter('vlambda', VLAMBDA)

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
        simulation.step(5000) #1ns

################################################################################
                ### NPT equilibration ###

        system.addForce(MonteCarloBarostat(1 * bar, 310 * kelvin))
        # to update the simulation object to take in account the new system
        simulation.context.reinitialize(True)
        print('Running NPT equilibration...')
        simulation.step(5000) #1ns

        # save the equilibration results to file
        simulation.saveState('equi.state')
        simulation.saveCheckpoint('equi.chk')

################################################################################
                ### Production run ###

        # Collect data
        xtc_reporter = XTCReporter('prod.xtc', 1000)
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
