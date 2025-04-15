import pytraj as pt
import argparse
import os

# Parse arguments
parser = argparse.ArgumentParser(description="Wrap and align trajectory using pytraj (method 2: no box).")
parser.add_argument("--dcd", required=True, help="Input DCD trajectory")
parser.add_argument("--pdb", required=True, help="Input PDB reference")
parser.add_argument("--stride", type=int, default=10, help="Stride interval for saving trajectory (default: 10)")
args = parser.parse_args()

# Load trajectory and reference
traj = pt.load(args.dcd, top=args.pdb)
ref = pt.load(args.pdb)

# Wrap with autoimage and align to reference
pt.autoimage(traj)
pt.align(traj, ref=ref[0], mask='@CA,C,N')

# Slice trajectory based on stride
traj_strided = traj[::args.stride]

# Output filenames
out_dcd = os.path.splitext(args.dcd)[0] + "_aligned.dcd"
out_pdb = os.path.splitext(args.dcd)[0] + "_aligned.pdb"

# Save strided aligned trajectory
pt.save(out_dcd, traj_strided, overwrite=True)
# Save first frame only
pt.save(out_pdb, traj_strided[:1], options="nobox", overwrite=True)

print(f"✅ Saved aligned trajectory (stride={args.stride}): {out_dcd}")
print(f"✅ Saved reference frame (no box): {out_pdb}")
