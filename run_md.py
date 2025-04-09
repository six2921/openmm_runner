import argparse
import os
import sys
import inspect
import json
import shutil

from datetime import datetime
from openmm.app import PDBFile, Simulation
from openmm import unit, XmlSerializer, LangevinIntegrator
from tqdm import tqdm
from openmm.app import StateDataReporter, DCDReporter
from pprint import pprint
from helpers import parse_time_data, extract_stats
from constants import integrator_config, forcefield_config, system_config
from reporter import items_energy, items_time

# --------------------------
# Argument Parsing
# --------------------------

parser = argparse.ArgumentParser(description="Run or append production MD simulation from checkpointed system.")
parser.add_argument("mode", choices=["run", "append"], help="Mode to start a new run or append to existing one")
parser.add_argument("-p", "--prev", required=True, help="Prefix of input files: expects .chk, .pdb, .xml (and .log/.time/.dcd/.json for append)")
parser.add_argument("-o", "--output", required=True, help="Base name for output files")
parser.add_argument("-t", "--time", type=float, required=True, help="Simulation time in nanoseconds")
parser.add_argument("--overwrite", action="store_true", help="Allow overwriting output files if they exist")
args = parser.parse_args()

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

# --------------------------
# File names to load and save
# --------------------------

input_prefix = args.prev
pdb_file = f"{input_prefix}.pdb"
xml_file = f"{input_prefix}.xml"
chk_file = f"{input_prefix}.chk"

output_prefix = args.output
json_output_file = f"{output_prefix}.json"
dcd_file = f"{output_prefix}.dcd"
log_file = f"{output_prefix}.log"
time_file = f"{output_prefix}.time"
out_pdb = f"{output_prefix}.pdb"
out_xml = f"{output_prefix}.xml"
out_chk = f"{output_prefix}.chk"

# --------------------------
# Validate Required Input Files & Prevent Overwrites
# --------------------------

required_files = [pdb_file, xml_file, chk_file]
if args.mode == "append":
    required_files += [
        f"{input_prefix}.log",
        f"{input_prefix}.time",
        f"{input_prefix}.dcd",
        f"{input_prefix}.json"
    ]

for file in required_files:
    if not os.path.exists(file):
        raise FileNotFoundError(f"Missing input file: {file}")

# Unified overwrite check (run & append)
output_files = [dcd_file, log_file, time_file, json_output_file]
if not args.overwrite:
    for file in output_files:
        if os.path.exists(file):
            raise FileExistsError(f"[ERROR] Output file '{file}' already exists. Use --overwrite to overwrite.")

# Append mode: copy previous outputs if input and output differ
if args.mode == "append" and input_prefix != output_prefix:
    shutil.copy(f"{input_prefix}.dcd", dcd_file)
    shutil.copy(f"{input_prefix}.log", log_file)
    shutil.copy(f"{input_prefix}.time", time_file)

# --------------------------
# Load PDB, Force Field System (.xml), and Checkpoint (.chk)
# --------------------------

pdb = PDBFile(pdb_file)
with open(xml_file) as f:
    system = XmlSerializer.deserialize(f.read())
with open(chk_file, "rb") as f:
    chk_data = f.read()

# --------------------------
# Set steps_interval and total_steps
# --------------------------

dt_ps = dt.value_in_unit(unit.picoseconds)

if args.mode == "run":
    report_interval_ps = 10 if args.time < 100 else 50
    steps_interval = int(report_interval_ps / dt_ps)
elif args.mode == "append":
    with open(f"{input_prefix}.json") as f:
        config_data = json.load(f)
    steps_interval = config_data["steps_interval"]
    report_interval_ps = steps_interval * dt_ps

total_steps = int(args.time * 1000 / dt_ps // 100 * 100)
chunk_steps = total_steps // 100

# --------------------------
# Initialize Simulation
# --------------------------

integrator = LangevinIntegrator(temp, gamma, dt)
simulation = Simulation(pdb.topology, system, integrator)

if args.mode == "run":
    dummy_sim = Simulation(pdb.topology, system, LangevinIntegrator(temp, gamma, dt))
    dummy_sim.context.setPositions(pdb.positions)
    dummy_sim.loadCheckpoint(chk_file)

    state = dummy_sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    positions = state.getPositions()
    velocities = state.getVelocities()
    box_vectors = state.getPeriodicBoxVectors()

    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)
    simulation.context.setPeriodicBoxVectors(*box_vectors)

elif args.mode == "append":
    simulation.context.setPositions(pdb.positions)
    simulation.loadCheckpoint(chk_file)

print(f"[INFO] Run Time : {args.time} ns ({total_steps} steps)")
print(f"[INFO] Trajectory Invertal : {report_interval_ps} ps")

# --------------------------
# StateDataReporter: .log(energy) and .time(time)
# --------------------------

simulation.reporters.append(DCDReporter(dcd_file, steps_interval, append=(args.mode == "append")))
simulation.reporters.append(StateDataReporter(
    time_file, steps_interval, append=(args.mode == "append"),
    **items_time, totalSteps=total_steps
))
simulation.reporters.append(StateDataReporter(
    log_file, steps_interval, append=(args.mode == "append"),
    **items_energy
))

# --------------------------
# Write .config file containing run configuration(start time, end time, etc.)
# --------------------------

# 1. Define shared config regardless of mode
common_config = {
    "output_prefix": output_prefix,
    "time_ns": args.time,
    "dt_ps": dt_ps,
    "total_steps": total_steps,
    "steps_interval": steps_interval,
    "steps_interval_ps": report_interval_ps,
    "started": datetime.now().strftime('%Y-%m-%d %H:%M:%S')
}

# 2. If append mode, calculate run extension values and expand config
append_config = {}
if args.mode == "append":
    time_segments, warnings = parse_time_data(time_file)
    latest = extract_stats(time_segments[max(time_segments)])

    start_step = int(latest.final_step)
    start_ns = float(start_step * dt_ps)
    end_step = start_step + total_steps
    end_ns = start_ns + args.time

    append_config = {
        "start_ns": start_ns,
        "end_ns": end_ns,
        "start_step": start_step,
        "end_step": end_step
    }

# 3. Combine and write final config
config_data = {
    **common_config,
    **append_config
}
with open(json_output_file, "w") as f:
    json.dump(config_data, f, indent=2)

# --------------------------
# Run Simulation
# --------------------------

print(f"[RUNNING] Running {args.time} ns production simulation......")
for _ in tqdm(range(100), desc="          Production"):
    simulation.step(chunk_steps)

# --------------------------
# Write outputs
# --------------------------

final_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)

with open(out_pdb, "w") as f:
    PDBFile.writeFile(simulation.topology, final_state.getPositions(), f)

with open(out_xml, "w") as f:
    f.write(XmlSerializer.serialize(simulation.system))

simulation.saveCheckpoint(out_chk)

print("\n[INFO] Production complete. Check outputs:")
print(f"  - {out_pdb}       # Final coordinates")
print(f"  - {dcd_file}       # Trajectory")
print(f"  - {log_file}       # Energy log")
print(f"  - {time_file}      # Progress log")
print(f"  - {out_chk}       # Checkpoint")
print(f"  - {json_output_file}      # Config")

