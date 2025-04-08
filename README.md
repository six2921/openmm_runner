## ⚠️ Caution
**This workflow does not sanitize protein or ligand structures.**
Make sure your input files are valid and clean before running simulations.

---

# 🧪 openmm_runner

A lightweight, modular Python-based workflow for running molecular dynamics (MD) simulations with [OpenMM](https://openmm.org/). Easily simulate apo or protein-ligand systems, monitor simulation progress, and align trajectories for visualization.

---

## 📁 Directory Structure and Script Roles

```
openmm_runner/
├── prepare_md.py         # System preparation (solvate, minimize, equilibrate)
├── run_md.py             # Production MD (new run or append mode)
├── check_time.py         # Check MD progress, speed, ETA
├── align_md.py           # Align and center trajectory using PyTraj
├── reporter.py           # Reporter config for energy/time logs
├── helpers.py            # Log/time parsing utilities
├── constants.py          # Default simulation parameters
├── system.pdb            # Final equilibrated structure (example)
```

### ▶️ Main Execution Scripts
- `prepare_md.py`: run solvation, minimization, and equilibration
- `run_md.py`: run or append production MD simulations

### ⚖️ Utility Tools
- `align_md.py`: align trajectory to a reference PDB using PyTraj
- `check_time.py`: parse `.time`/`.json` files and show progress, speed, ETA

### 🎓 Support Modules
- `reporter.py`: preconfigured reporter items for OpenMM
- `constants.py`: force field and integrator settings
- `helpers.py`: utility functions for logging, parsing

---

## 🚀 How to Use

### 1️⃣ Prepare System

#### Protein-only (apo)
```bash
python prepare_md.py -p protein.pdb
```

#### Protein-ligand complex
```bash
python prepare_md.py -p protein.pdb -l ligand.sdf
```

> Ligand file must be `.mol` or `.sdf`.

**Output files:**
```
system_initial.pdb      # unsolvated system
system_solvated.pdb     # after solvation
system_minimized.pdb    # after energy minimization
system_nvt.log          # NVT equilibration
system_npt.log          # NPT equilibration
system.pdb              # final coordinates
system.xml              # serialized OpenMM system
system.chk              # checkpoint
```

---

### 2️⃣ Run Production MD
```bash
python run_md.py run -p system -o md -t 100
```
- `-t` is simulation time in **nanoseconds**

**Output files:**
```
md.dcd              # trajectory
md.pdb              # final coordinates
md.chk              # checkpoint
md.log              # energy log
md.time             # progress log
md.json             # run configuration
```

---

### 3️⃣ Append Simulation
```bash
python run_md.py append -p md -o md2 -t 50
```
- `-t` is in **nanoseconds**
- `-p` and `-o` **must not be the same**
- The script copies all files from `-p` and appends results to the new set
- If you want to reuse the old name, **delete the original `.dcd` manually** first

---

### 4️⃣ Align Trajectory
```bash
python align_md.py --dcd md.dcd --pdb system_solvated.pdb
```
- Wraps and aligns trajectory using protein backbone atoms (`@CA,C,N`)
- Overwrites the original `.dcd` and `.pdb` files

---

### 5️⃣ Check Simulation Progress
```bash
python check_time.py md
```
Sample output:
```
Progress (ns)       : 0.1 / 0.1 (100%)
Steps Range (Steps) : 200000 ~ 250000
Time Range (ns)     : 0.4 ~ 0.5
Speed (ns/day)      : 562.1
Estimated Duration  : 0d 0h 0m 15s
Estimated End       : 2025-04-08 14:56:37
```

---

## ✨ Features
- Run full OpenMM MD simulation with only protein and ligand files, no complex setup
- Automatically parameterize ligands using **GAFF**
- Supports **append mode** to continue previous runs
- Includes tool to **check simulation progress, speed, ETA**

---

## 🔗 License
MIT License. Modify and use freely. Not intended for commercial redistribution.


