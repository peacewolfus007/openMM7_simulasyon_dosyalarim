from openmm import unit
import numpy as np
import openmm as mm
from openmm import app

# Simulation parameters
n_iterations = 250
n_steps_per_iteration = 250
exchange_attempt_frequency = 100

# Temperature setup
temperatures = [300*unit.kelvin, 310*unit.kelvin, 320*unit.kelvin, 330*unit.kelvin]

# Boltzmann constant in kJ/(mol*K)
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA

# Setup simulations
def setup_simulation(temp, prmtop_file, crd_file):
    # Load the Amber files
    prmtop = app.AmberPrmtopFile(prmtop_file)
    inpcrd = app.AmberInpcrdFile(crd_file)
    
    # Create the system
    system = prmtop.createSystem(nonbondedMethod=app.PME, 
                                 nonbondedCutoff=1*unit.nanometer,
                                 constraints=app.HBonds)
    
    # Create an integrator
    integrator = mm.LangevinIntegrator(temp, 1/unit.picosecond, 0.002*unit.picoseconds)
    
    # Create the simulation
    simulation = app.Simulation(prmtop.topology, system, integrator)
    
    # Set initial positions
    simulation.context.setPositions(inpcrd.positions)
    
    # Check if velocities are present and set them
    if inpcrd.velocities is not None:
        simulation.context.setVelocities(inpcrd.velocities)
    else:
        # If no velocities in inpcrd, generate them
        simulation.context.setVelocitiesToTemperature(temp)
    
    return simulation

# Amber file paths
prmtop_file = 'fizo.prmtop'
crd_file = 'fizo.crd'

# Create simulations for each temperature
simulations = [setup_simulation(temp, prmtop_file, crd_file) for temp in temperatures]

# Simulation loop
for iteration in range(n_iterations):
    for simulation in simulations:
        simulation.step(n_steps_per_iteration)

    if iteration % exchange_attempt_frequency == 0:
        for i in range(len(temperatures) - 1):
            state1 = simulations[i].context.getState(getEnergy=True)
            state2 = simulations[i+1].context.getState(getEnergy=True)

            E1 = state1.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            E2 = state2.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

            T1 = temperatures[i].value_in_unit(unit.kelvin)
            T2 = temperatures[i+1].value_in_unit(unit.kelvin)

            delta = (E2 - E1) * (1/(kB * T1) - 1/(kB * T2))

            if delta < 0 or np.random.rand() < np.exp(-delta):
                positions1, velocities1 = state1.getPositions(), state1.getVelocities()
                positions2, velocities2 = state2.getPositions(), state2.getVelocities()

                simulations[i].context.setPositions(positions2)
                simulations[i].context.setVelocities(velocities2)
                simulations[i+1].context.setPositions(positions1)
                simulations[i+1].context.setVelocities(velocities1)

print("Simulations completed.")