import openmm as mm
from openmm import app
from openmm.unit import *
from mdtraj.reporters import NetCDFReporter
from openmm.app import StateDataReporter 
from sys import stdout
import time

# Dosyaları yükleme
prmtop = app.AmberPrmtopFile('fizo.prmtop')
inpcrd = app.AmberInpcrdFile('fizo.crd')

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
for step in range(5000): #(her adım 2 fs)
    simulation.step(2)

# Veri kaydetme
simulation.reporters.append(NetCDFReporter('fizo.nc', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# MD üretimi
print('Running Production...')
start_time = time.time()  # Başlangıç zamanını kaydet
last_checkpoint_time = start_time
num_steps_per_ns = 500000  # 2 fs per step olduğundan, 1 ns için 500,000 adım
checkpoint_interval = 2 * num_steps_per_ns  # Her 2 ns için

for step in range(0, 1000000, checkpoint_interval):
    simulation.step(checkpoint_interval)
    current_time = time.time()
    elapsed_time = current_time - last_checkpoint_time
    last_checkpoint_time = current_time
    checkpoint_file = f'fizo_{step + checkpoint_interval}.chk'
    simulation.saveCheckpoint(checkpoint_file)
    print(f'Checkpoint saved to {checkpoint_file}, Time since last checkpoint: {elapsed_time:.2f} seconds')

print('Done!')