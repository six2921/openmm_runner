# 🧬 openmm_runner

A clean, modular Python-based workflow for running molecular dynamics (MD) simulations using [OpenMM](https://openmm.org/).  
Includes everything from system preparation to production runs, trajectory alignment, and progress tracking.

---

## 📦 Repository Structure

```
openmm_runner/
├── prepare_md.py         # System preparation (solvate, minimize, equilibrate)
├── run_md.py             # Production MD (new run or append mode)
├── align_md.py           # Align and center trajectory using PyTraj
├── check_time.py         # Check MD progress, speed, ETA
├── reporter.py           # Reporter config for energy/time logs
├── helpers.py            # Log/time parsing utilities
├── constants.py          # Default simulation parameters
├── system.pdb            # Final equilibrated structure (example)
├── md.dcd                # Example trajectory
├── md.pdb                # Final snapshot
├── md.log                # Energy log
├── md.time               # Progress log
├── md.json               # Run config metadata
```

---

## 🔧 Requirements

Install all dependencies with Conda:

```bash
conda create -n openmm python=3.10 -c conda-forge openmm pytraj tqdm numpy
conda activate openmm
```

---

## 🧪 How to Use

### 1️⃣ Prepare System

#### Protein-only (apo)

```bash
python prepare_md.py -p protein.pdb
```

#### Protein-ligand complex

```bash
python prepare_md.py -p protein.pdb -l ligand.sdf
```

📂 Output files:

```
  - system_initial.pdb
  - system_solvated.pdb
  - system_minimized.pdb
  - system_nvt.log
  - system_npt.log
  - system.pdb
  - system.xml
  - system.chk
```

---

### 2️⃣ Run Production MD

```bash
python run_md.py run -p system -o md -t 100
```

- `-t` is in nanoseconds
- Auto-aligns trajectory to `system_solvated.pdb` by default

📂 Output files:

```
  - md.pdb           # Final coordinates
  - md.dcd           # Trajectory
  - md.log           # Energy log
  - md.time          # Progress log
  - md.chk           # Checkpoint
  - md.json          # Config (for append)
```

---

### 3️⃣ Append Simulation

```bash
python run_md.py append -p md -o md_extension -t 100
```

- Appends 100 ns to `md` simulation
- Still aligns output trajectory unless disabled

---

### 4️⃣ Align Trajectory Manually

```bash
python align_md.py --dcd md.dcd --pdb system_solvated.pdb
```

- Aligns trajectory to the reference
- Autoimages and centers the protein
- Overwrites original `.dcd` and `.pdb`

---

### 5️⃣ Check Simulation Progress

```bash
python check_time.py md
```

🧾 Output example:

```
Progress (ns)       : 0.1 / 0.1 (100%)
Steps Range (Steps) : 200000 ~ 250000
Time Range (ns)     : 0.4 ~ 0.5
Speed (ns/day)      : 562.1
Estimated Duration  : 0d 0h 0m 15s
Estimated End       : 2025-04-08 14:56:37
```

---

## ✅ Features

- Protein-only or protein-ligand system setup
- Appendable production MD with checkpointing
- PyTraj-based alignment and centering
- Config logging and progress tracking
- Safe overwrite protection and clean interface

---

## 🔗 License

MIT License.  
Use, modify, and extend freely for your own MD workflows and research.

---

## 💬 Contact

For questions or suggestions, feel free to open an issue or pull request.
