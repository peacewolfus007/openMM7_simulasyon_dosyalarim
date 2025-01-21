import openmm as mm
from openmm import app
from openmm.unit import *
from mdtraj.reporters import NetCDFReporter
from openmm.app import StateDataReporter 
from sys import stdout

# Dosyaları yükleme
prmtop = app.AmberPrmtopFile('kaempfexx.prmtop')
inpcrd = app.AmberInpcrdFile('kaempfexx.crd')

# Sistemi oluşturma
system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometers, constraints=app.HBonds)
integrator = mm.LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = app.Simulation(prmtop.topology, system, integrator)

# Koordinatları ayarlama
simulation.context.setPositions(inpcrd.positions)

# Enerji minimizasyonu
print('Minimizing...')
simulation.minimizeEnergy(maxIterations=500)

# Isıtma
print('Heating...')
simulation.context.setVelocitiesToTemperature(300*kelvin)
for step in range(500000): #(her adım 2 fs)
    simulation.step(2)

# Veri kaydetme
simulation.reporters.append(NetCDFReporter('kaempfexx.nc', 1000))
simulation.saveCheckpoint('kaempfexx.chk')
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# MD üretimi
print('Running Production...')
simulation.step(1000000) # 0.1 ns için 50,000 adım
print('Done!')
