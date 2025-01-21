import openmm as mm
from openmm import app
from openmm.unit import *
from mdtraj.reporters import NetCDFReporter
from openmm.app import StateDataReporter 
from sys import stdout
import time
import shutil
import os

# Dosyaları yükleme
prmtop = app.AmberPrmtopFile('fizo_hmass.prmtop')
inpcrd = app.AmberInpcrdFile('fizo.crd')

# Sistemi oluşturma
system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometers, constraints=app.HBonds)
integrator = mm.LangevinIntegrator(300*kelvin, 1/picosecond, 4*femtoseconds)
simulation = app.Simulation(prmtop.topology, system, integrator)

# Koordinatları ayarlama
simulation.context.setPositions(inpcrd.positions)

# Enerji minimizasyonu
print('Minimizing...')
simulation.minimizeEnergy(maxIterations=500)

# Isıtma
print('Heating...')
simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.step(10000)  # 10000 adım = 20 ps

# Veri kaydetme
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# MD üretimi
print('Running Production...')
total_steps = 50000  # 100 ns için toplam adım sayısı (2 fs/adım)
checkpoint_interval = 1000  # Her 2 ns için adım sayısı

for step in range(0, total_steps, checkpoint_interval):
    start_time = time.time()
    
    current_ns = (step + checkpoint_interval) * 4 / 1000  # ns cinsinden mevcut zaman
    nc_file = f'fizo_{current_ns:.0f}ns.nc'
    chk_file = f'fizo_{current_ns:.0f}ns.chk'
    
    # Her döngüde yeni bir NetCDF reporter ekleyin
    nc_reporter = NetCDFReporter(nc_file, checkpoint_interval)
    simulation.reporters.append(nc_reporter)
    
    simulation.step(checkpoint_interval)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    simulation.saveCheckpoint(chk_file)
    
    print(f'{current_ns:.0f} ns tamamlandı. Süre: {elapsed_time:.2f} saniye')
    print(f'NetCDF dosyası kaydedildi: {nc_file}')
    print(f'Checkpoint dosyası kaydedildi: {chk_file}')
    
    # NetCDF reporter'ı kaldırın
    simulation.reporters.pop()
    
    # Drive'a kaydet (eğer Google Colab kullanıyorsanız)
    drive_path = '/content/drive/My Drive/md_simulation/'
    if os.path.exists(drive_path):
        shutil.copy(nc_file, drive_path)
        shutil.copy(chk_file, drive_path)
        print(f'Dosyalar Drive\'a kopyalandı: {nc_file}, {chk_file}')

print('Simülasyon tamamlandı!')