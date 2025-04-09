### ⚠️ Caution
**This workflow does not sanitize protein or ligand structures.**
Make sure your input files are valid and clean before running simulations.

---

### 🧪 openmm_runner

A lightweight, modular Python-based workflow for running molecular dynamics (MD) simulations with [OpenMM](https://openmm.org/). Easily run full OpenMM MD simulations for protein-ligand systems in just two commands.

📦 Installation
```bash
conda create -n openmm python=3.10
conda activate openmm
conda install -c conda-forge openmm openff-toolkit rdkit mdtraj tqdm
```

---

### 📁 Scripts

**MD Runner**
- `prepare_md.py`: run solvation, minimization, and equilibration
- `run_md.py`: run or append production MD simulations

**Align&Monitoring**
- `align_md.py`: align trajectory to a reference PDB using PyTraj
- `check_time.py`: parse `.time`/`.json` files and show progress, speed, ETA

**Configurations**
- `reporter.py`: preconfigured reporter items for OpenMM
- `constants.py`: force field and integrator settings
- `helpers.py`: utility functions for logging, parsing

---

### 🚀 How to Use

**1️⃣ Prepare System**

**Protein-only**
```bash
python prepare_md.py -p protein.pdb
```

**Protein-ligand complex**
```bash
python prepare_md.py -p protein.pdb -l ligand.sdf
```
- Ligand file must be `.mol` or `.sdf`.

**Output files:**
```
system_initial.pdb      # system before solvation
system_solvated.pdb     # after solvation
system_minimized.pdb    # after energy minimization
system_nvt.log          # NVT equilibration
system_npt.log          # NPT equilibration
system_relax.pdb        # final coordinates
system_relax.xml        # serialized OpenMM system
system_relax.chk        # checkpoint
```

**2️⃣ Run Production MD**
```bash
python run_md.py run -p system_relax -o md -t 100
```
- `-p` is basename of pdb, xml, and chk files of system preparation
- `-t` is simulation time in **nanoseconds**

**Output files:**
```
md.dcd              # trajectory
md.pdb              # final coordinates
md.chk              # checkpoint
md.log              # energy log
md.time             # progress log
md.json             # run information
```

**3️⃣ Check Simulation Progress**
```bash
python check_time.py md
```
- Specify the base name of the production run to check its progress

Sample output:
```
Progress (ns)       : 0.1 / 0.1 (100%)
Steps Range (Steps) : 200000 ~ 250000
Time Range (ns)     : 0.4 ~ 0.5
Speed (ns/day)      : 562.1
Estimated Duration  : 0d 0h 0m 15s
Estimated End       : 2025-04-08 14:56:37
```

**4️⃣ Align Trajectory**
```bash
python align_md.py --dcd md.dcd --pdb system_solvated.pdb
```
- --dcd is trajectory you want to align. _aligned.dcd will be created.
- --pdb is reference protein for alignment. system_solvated.pdb is recommended.

**5️⃣ Append Simulation**
```bash
python run_md.py append -p md -o md2 -t 50
```
- -p means previous. Input basename of the previous run. eg. md
- `-o` 	must not be the same as an existing file name. Use --overwrite to allow this.

---

## 🔗 License
MIT License. Modify and use freely. Not intended for commercial redistribution.
