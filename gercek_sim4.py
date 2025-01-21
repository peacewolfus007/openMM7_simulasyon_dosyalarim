import openmm as mm
from openmm import app
from openmm.unit import *
from mdtraj.reporters import NetCDFReporter
from sys import stdout
import mdtraj as md

# İlk simülasyondan son durumu yükle
trajectory = md.load('output3.nc', top='complex_su.prmtop')
last_frame = trajectory[-1]

# Dosyaları yükleme
prmtop = app.AmberPrmtopFile('complex_su.prmtop')
inpcrd = app.AmberInpcrdFile('complex_su.crd')

# Sistemi oluşturma
system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometers, constraints=app.HBonds)
integrator = mm.LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = app.Simulation(prmtop.topology, system, integrator)
# İlk simulasyondan hızları yükle
simulation.loadCheckpoint('checkpoint3.chk')
simulation.saveCheckpoint('checkpoint4.chk')

# Veri kaydetme
simulation.reporters.append(NetCDFReporter('output4.nc', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# MD üretimi
print('Running Production...')
simulation.step(500000) # 0.1 ns için 50,000 adım
print('Done!')
