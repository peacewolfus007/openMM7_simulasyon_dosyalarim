import openmm as mm
from openmm import app
from openmm.unit import *
from mdtraj.reporters import NetCDFReporter
from sys import stdout
import mdtraj as md

# Dosyaları yükleme
prmtop = app.AmberPrmtopFile('complex_su.prmtop')

# Sistemi oluşturma
system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometers, constraints=app.HBonds)
integrator = mm.LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = app.Simulation(prmtop.topology, system, integrator)

# En son checkpoint dosyasını yükleme
simulation.loadCheckpoint('checkpoint2.chk')

# Veri kaydetme
simulation.reporters.append(NetCDFReporter('output4.nc', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# MD üretimi
print('Running Production...')
simulation.step(500000) # 1 ns için 500,000 adım
print('Done!')
