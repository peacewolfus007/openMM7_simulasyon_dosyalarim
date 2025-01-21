import openmm as mm
from openmm import app
from openmm.unit import *
from openmm.app import StateDataReporter, DCDReporter
from sys import stdout

# Dosyaları yükleme
prmtop = app.AmberPrmtopFile('fizo.prmtop')
inpcrd = app.AmberInpcrdFile('fizo.crd')

# Sistemi oluşturma
system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometers, constraints=app.HBonds)
integrator = mm.LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtoseconds)  # Zaman adımını 1 fs olarak ayarla
simulation = app.Simulation(prmtop.topology, system, integrator)

# Koordinatları ayarlama
simulation.context.setPositions(inpcrd.positions)

# Enerji minimizasyonu
print('Minimizing...')
simulation.minimizeEnergy(maxIterations=350)

# Isıtma
print('Heating...')
initial_temp = 100*kelvin
temp_increment = 10*kelvin
while initial_temp < 300*kelvin:
    simulation.context.setVelocitiesToTemperature(initial_temp)
    initial_temp += temp_increment
    simulation.step(100)  # Her adım 1 fs

# Veri kaydetme
simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.saveCheckpoint('fizo.chk')
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# MD üretimi
print('Running Production...')
simulation.step(50000)  # 0.1 ns için 50,000 adım
print('Done!')
