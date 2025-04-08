#!/usr/bin/env python
import argparse
import csv
from datetime import datetime, timedelta
import os
import json


def format_duration(seconds):
    """Convert seconds to 'Xd Xh Xm Xs' string."""
    seconds = int(round(seconds))
    d, rem = divmod(seconds, 86400)
    h, rem = divmod(rem, 3600)
    m, s = divmod(rem, 60)
    return f"{d}d {h}h {m}m {s}s"


def format_datetime(dt):
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def format_ns(value):
    return str(int(round(value))) if value >= 10 else f"{value:.4g}"


# --- Argument Parsing ---
parser = argparse.ArgumentParser(
    description="Check MD simulation progress, speed, and estimated finish time."
)
parser.add_argument("basename", help="Base name for .json and .time files (e.g., md4)")
args = parser.parse_args()

# --- File Paths ---
base = args.basename
json_file = f"{base}.json"
time_file = f"{base}.time"

if not os.path.exists(json_file):
    raise FileNotFoundError(f"Config file not found: {json_file}")
if not os.path.exists(time_file):
    raise FileNotFoundError(f"Time file not found: {time_file}")

# --- Read Config File (.json) ---
with open(json_file, "r") as f:
    config = json.load(f)

try:
    total_steps = int(config["total_steps"])
    dt_ps = float(config["dt_ps"])              # step size in ps
    target_ns = float(config["time_ns"])        # total MD for this run
    steps_per_interval = int(config["steps_interval"])
    started_str = config["started"]
    start_step = int(config.get("start_step", 0))  # default to 0 if not present
    end_step = int(config.get("end_step", total_steps))
except KeyError as e:
    raise KeyError(f"Missing key in config file: {e}")

overall_start = datetime.strptime(started_str, "%Y-%m-%d %H:%M:%S")

# --- Read Time File ---
with open(time_file, "r") as f:
    reader = csv.DictReader(f)
    records = list(reader)

if not records:
    raise RuntimeError("Time file is empty.")

# --- Identify start and end steps for this MD run ---
last_dash_index = max(
    i for i, row in enumerate(records)
    if row.get("Time Remaining", "").strip() == "--"
)
start_record = records[last_dash_index]
start_elapsed = float(start_record["Elapsed Time (s)"])

# --- Latest record = current step ---
end_record = records[-1]
current_step = int(end_record["Step"])
current_elapsed = float(end_record["Elapsed Time (s)"])

# --- Compute progress ---
steps_done = current_step - start_step
done_ns = steps_done * dt_ps / 1000
total_ns = total_steps * dt_ps / 1000
steps_left = total_steps - steps_done

time_spent = current_elapsed - start_elapsed
if steps_done > 0:
    time_per_step = time_spent / steps_done
    estimated_total_time = total_steps * time_per_step
    remaining_time = estimated_total_time - time_spent
    ns_per_day = done_ns * 86400 / time_spent
else:
    time_per_step = estimated_total_time = remaining_time = ns_per_day = 0

# --- Time Stamps ---
current_actual_time = overall_start + timedelta(seconds=current_elapsed)
estimated_end_time = current_actual_time + timedelta(seconds=remaining_time)

# --- Range in ns ---
range_start_step = start_step
range_end_step = start_step + total_steps
range_start_ns = range_start_step * dt_ps / 1000
range_end_ns = range_end_step * dt_ps / 1000

# --- Final Output ---
print(f" ")
print(f"Progress (ns)       : {format_ns(done_ns)} / {format_ns(total_ns)} ({done_ns / total_ns * 100:.0f}%)")
print(f"Steps Range (Steps) : {range_start_step:,} ~ {range_end_step:,}")
print(f"Time Range (ns)     : {format_ns(range_start_ns)} ~ {format_ns(range_end_ns)}")
print(f"Speed (ns/day)      : {ns_per_day:.1f}")
print(f"Estimated Duration  : {format_duration(estimated_total_time)}")
print(f"Estimated End       : {format_datetime(estimated_end_time)}")
print(f" ")
