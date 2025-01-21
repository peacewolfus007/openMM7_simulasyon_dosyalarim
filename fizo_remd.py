import openmm as mm
from openmm import app
from openmm.unit import *
from sys import stdout, setrecursionlimit
import numpy as np

# Özyineleme derinliğini artırma
setrecursionlimit(1500)

# Dosyaları yükleme
prmtop = app.AmberPrmtopFile('fizo.prmtop')
inpcrd = app.AmberInpcrdFile('fizo.crd')

# Sıcaklıkları belirleme
temperatures = [269.5, 300.0, 334.0, 371.8, 413.9, 460.7, 512.9, 570.9]*kelvin
n_replicas = len(temperatures)

# Sistem ve entegratör durumlarını oluşturma
systems = [prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometers, constraints=app.HBonds) for _ in temperatures]
integrators = [mm.LangevinIntegrator(temp, 1/picosecond, 1*femtosecond) for temp in temperatures]

# Simülasyonları oluşturma
simulations = [app.Simulation(prmtop.topology, system, integrator) for system, integrator in zip(systems, integrators)]

def check_positions_for_nan(positions):
    for pos in positions:
        if np.isnan(pos.value_in_unit(nanometers)).any():
            return True
    return False

# Koordinatları ayarlama ve enerji minimizasyonu
for simulation in simulations:
    simulation.context.setPositions(inpcrd.positions)
    simulation.minimizeEnergy(maxIterations=5000)  # Maksimum iterasyonu artırın
    positions = simulation.context.getState(getPositions=True).getPositions()
    if check_positions_for_nan(positions):
        raise ValueError("NaN detected in minimized positions.")
    print(f"Minimized positions for replica: {positions}")

# Simülasyonları başlatma
n_steps_per_iteration = 250
n_iterations = 250
exchange_attempt_frequency = 100  # her 100 adımda bir değişim denemesi
print(f"Simülasyon başlatılıyor: {n_iterations} iterasyon, her iterasyonda {n_steps_per_iteration} adım.")

kB = BOLTZMANN_CONSTANT_kB  # Boltzmann sabitini uygun birime dönüştür

for iteration in range(n_iterations):
    for simulation in simulations:
        simulation.step(n_steps_per_iteration)
        positions = simulation.context.getState(getPositions=True).getPositions()
        if check_positions_for_nan(positions):
            raise ValueError("NaN detected in particle positions during simulation.")
    
    # Enerji hesaplama ve replika değişim denemesi
    for i in range(n_replicas - 1):
        state1 = simulations[i].context.getState(getEnergy=True)
        state2 = simulations[i+1].context.getState(getEnergy=True)
        E1 = state1.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        E2 = state2.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        T1 = temperatures[i].value_in_unit(kelvin)
        T2 = temperatures[i+1].value_in_unit(kelvin)
        
        delta = (E2 - E1) * (1/(kB * T1) - 1/(kB * T2))
        delta_value = delta  # Birimleri zaten uyumlu, doğrudan delta'yı kullanın
        
        if delta_value < 0 or np.random.rand() < np.exp(-delta_value):
            # Replikaların pozisyonlarını ve hızlarını değiştir
            positions1 = state1.getPositions(asNumpy=True)
            velocities1 = state1.getVelocities(asNumpy=True)
            positions2 = state2.getPositions(asNumpy=True)
            velocities2 = state2.getVelocities(asNumpy=True)
            
            simulations[i].context.setPositions(positions2)
            simulations[i].context.setVelocities(velocities2)
            simulations[i+1].context.setPositions(positions1)
            simulations[i+1].context.setVelocities(velocities1)
    print("Replika değişim denemesi tamamlandı.")

    # Checkpoint ve çıktı dosyaları
    for i, temp in enumerate(temperatures):
        simulations[i].saveCheckpoint(f'hesperetin_{int(temp._value)}K.chk')
        simulations[i].reporters.append(NetCDFReporter(f'trajectory_{int(temp._value)}K.nc', 1000))
        simulations[i].reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
print("Simülasyon tamamlandı.")
