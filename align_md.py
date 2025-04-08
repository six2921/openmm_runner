import pytraj as pt
import argparse
import os

parser = argparse.ArgumentParser(description="Wrap and align trajectory using pytraj (method 2: no box).")
parser.add_argument("--dcd", required=True)
parser.add_argument("--pdb", required=True)
args = parser.parse_args()

# Load trajectory and reference structure
traj = pt.load(args.dcd, top=args.pdb)
ref = pt.load(args.pdb)

# PBC wrapping: center protein in box
pt.autoimage(traj)

# Align to reference
pt.align(traj, ref=ref[0], mask='@CA,C,N')

# Save aligned trajectory and first frame
out_dcd = os.path.splitext(args.dcd)[0] + "_aligned.dcd"
out_pdb = os.path.splitext(args.dcd)[0] + "_aligned.pdb"

pt.save(out_dcd, traj, overwrite=True)
pt.save(out_pdb, traj[:1], options="nobox", overwrite=True)

print(f"✅ Saved aligned trajectory: {out_dcd}")
print(f"✅ Saved reference frame (no box): {out_pdb}")

