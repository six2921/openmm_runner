import argparse
import os
from openmm.app import PDBFile, Modeller, Simulation, StateDataReporter, HBonds
from openmm import unit, Platform, MonteCarloBarostat, LangevinIntegrator, XmlSerializer
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
from tqdm import tqdm
from pprint import pprint
from constants import integrator_config, forcefield_config, system_config
from reporter import items_energy, items_time

# --------------------------
# Argument Parsing
# --------------------------

parser = argparse.ArgumentParser(description="Prepare solvated, minimized, and equilibrated system for production MD")
parser.add_argument("-p", "--protein", required=True, help="Protein PDB file")
parser.add_argument("-l", "--ligand", required=False, help="Ligand file (.mol or .sdf)")
parser.add_argument("-o", "--output", default="system_relax", help="Base name for output files")
parser.add_argument("--overwrite", action="store_true", help="Allow overwriting output files if they exist")
args = parser.parse_args()

# --------------------------
# Overwrite protection
# --------------------------
if not args.overwrite:
    for ext in [".dcd", ".pdb", ".xml", ".chk"]:
        output_file = f"{args.output}{ext}"
        if os.path.exists(output_file):
            raise FileExistsError(f"[ERROR] Output file '{output_file}' already exists. Use --overwrite to overwrite.")

# ---------------
# Configuration unpacking
# ---------------

dt = integrator_config["dt"]
temp = integrator_config["temp"]
gamma = integrator_config["gamma"]

pressure = system_config["pressure"]
baro_temp = system_config["baro_temp"]
baro_freq = system_config["baro_freq"]

# Print configuration
print("  ---------- [ Integrator ] ----------")
pprint(integrator_config)
print("")

print("  ---------- [ System ] ----------")
pprint(system_config)
print("")

print("  ---------- [ Forcefield ] ----------")
pprint(forcefield_config)
print("")

# ---------------
# Input loading
# ---------------

protein = PDBFile(args.protein)
ligand = Molecule.from_file(args.ligand) if args.ligand else None

# Set molecule-related kwargs only if ligand is provided
small_molecule_args = {
    "small_molecule_forcefield": forcefield_config["small_molecule_ff"],
    "molecules": [ligand]
} if ligand else {}

system_generator = SystemGenerator(
    forcefields=forcefield_config["forcefields"],
    forcefield_kwargs=forcefield_config["forcefield_kwargs"],
    **small_molecule_args
)

# Topology generation
modeller = Modeller(protein.topology, protein.positions)
if ligand:
    modeller.add(
        ligand.to_topology().to_openmm(),
        ligand.conformers[0].to_openmm())

# Save initial system
with open("system_initial.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

# ---------------
# Solvation
# ---------------

print("[RUNNING] Solvation......")
modeller.addSolvent(
    system_generator.forcefield,
    padding=1.0 * unit.nanometer,
    ionicStrength=0.15 * unit.molar,
    neutralize=True)

# Save solvated structure
with open("system_solvated.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

# ---------------
# System setup
# ---------------

system = system_generator.create_system(
    modeller.topology,
    molecules=[ligand] if ligand else [])

integrator = LangevinIntegrator(temp, gamma, dt)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# ---------------
# Minimization
# ---------------

print("[RUNNING] Minimization......")
simulation.minimizeEnergy()
positions = simulation.context.getState(getPositions=True).getPositions()

with open("system_minimized.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, positions, f)

# ---------------
# NVT Equilibration
# ---------------

time_nvt = 100  # total time (ps)
print(f"[RUNNING] Running {time_nvt} ps NVT equilibration......")

simulation.reporters.append(StateDataReporter("system_nvt.log", 1000, **items_energy))
simulation.context.setVelocitiesToTemperature(temp)

for _ in tqdm(range(100), desc="          NVT"):
    simulation.step(int(time_nvt / (100 * dt.value_in_unit(unit.picoseconds))))

import time  # 파일 상단에서 import 필요
import datetime  # ETA 포맷용

# ---------------
# NPT Equilibration
# ---------------

time_npt = 100  # total time (ps)
print(f"[RUNNING] Running {time_npt} ps NPT equilibration......")

# ---Get current state from NVT ---
state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
positions = state.getPositions()
box_vectors = state.getPeriodicBoxVectors()

# --- Add barostat and re-create simulation ---
system.addForce(MonteCarloBarostat(pressure, baro_temp, baro_freq))
integrator = LangevinIntegrator(temp, gamma, dt)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(positions)
simulation.context.setPeriodicBoxVectors(*box_vectors)
simulation.context.setVelocitiesToTemperature(temp)

# --- Set reporter and run ---
simulation.reporters.append(StateDataReporter("system_npt.log", 1000, **items_energy))

# --- Timing ---
start_time = time.time()

for _ in tqdm(range(100), desc="          NPT"):
    simulation.step(int(time_npt / (100 * dt.value_in_unit(unit.picoseconds))))

elapsed = time.time() - start_time
estimated_100ns = elapsed * (100_000 / time_npt)  # 100,000 ps = 100 ns
eta = str(datetime.timedelta(seconds=int(estimated_100ns)))   # ETA will be printed at the last

# ---------------
# Final output
# ---------------

final_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)

with open(f"{args.output}.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, final_state.getPositions(), f)  # Final coordinates (after NPT)

with open(f"{args.output}.xml", "w") as f:
    f.write(XmlSerializer.serialize(simulation.system))  # Final system state (forces, constraints, barostat, etc.)

simulation.saveCheckpoint(f"{args.output}.chk")  # Checkpoint to resume simulation

print(" ")
print("[INFO] System prepared. Check outputs:")
print("  - system_initial.pdb")
print("  - system_solvated.pdb")
print("  - system_minimized.pdb")
print("  - system_nvt.log")
print("  - system_npt.log")
print(f"  - {args.output}.pdb")
print(f"  - {args.output}.xml")
print(f"  - {args.output}.chk")

print()
print("[INFO] Estimated time for 100 ns simulation (based on NPT performance):", eta)
