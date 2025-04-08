import os, sys, inspect
from datetime import datetime
from openmm.app import StateDataReporter, DCDReporter, PDBFile, Simulation
from openmm import LangevinIntegrator, unit, XmlSerializer

# --------------------------------------------------
# parse_time_data: Parse .time file
# --------------------------------------------------

from dataclasses import dataclass
import pandas as pd
from io import StringIO

@dataclass
class TimeStats:
    interval_steps: int
    final_step: int
    first_step: int
    total_steps: int
    total_elapsed: float
    time_per_step: float

def parse_time_data(file_path):
    """
    Parse the MD time log file and split the data by segments separated by "--".
    
    Returns:
        A tuple containing:
          - A dictionary mapping segment numbers to pandas DataFrames.
          - A list of warning messages (if any required columns are missing or intervals are inconsistent).
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # 헤더에서 주석 제거
    if lines[0].startswith("#"):
        lines[0] = lines[0][1:]
    content = ''.join(lines)
    df = pd.read_csv(StringIO(content))

    # 필수 컬럼 확인
    required_columns = ["Step", "Time (ps)", "Elapsed Time (s)", "Time Remaining"]
    missing_columns = [col for col in required_columns if col not in df.columns]
    warnings = []
    if missing_columns:
        warnings.append(f"Missing required columns: {missing_columns}")

    # Step 간격이 일정한지 확인
    step_diffs = df["Step"].diff().dropna()
    if not (step_diffs == step_diffs.iloc[0]).all():
        warnings.append("Warning: Step intervals are not consistent")

    # "--" 기준 데이터 분할
    split_indices = df.index[df["Time Remaining"] == "--"].tolist()
    split_datasets = {}
    if split_indices:
        prev_idx = 0
        count = 1
        for idx in split_indices + [len(df)]:
            part = df.iloc[prev_idx:idx+1].copy()
            part.columns.name = None
            split_datasets[count] = part
            prev_idx = idx + 1
            count += 1
    else:
        split_datasets[0] = df

    return split_datasets, warnings

def extract_stats(data: pd.DataFrame) -> TimeStats:
    """
    단일 세그먼트에 대해 통계값을 계산해 TimeStats 객체로 반환
    """
    interval = data["Step"].iloc[1] - data["Step"].iloc[0]
    final_step = data["Step"].iloc[-1]
    first_step = data["Step"].iloc[0] - interval
    total_steps = final_step - first_step
    total_elapsed = data["Elapsed Time (s)"].iloc[-1]
    time_per_step = total_elapsed / total_steps if total_steps != 0 else None

    return TimeStats(
        interval_steps=interval,
        final_step=final_step,
        first_step=first_step,
        total_steps=total_steps,
        total_elapsed=total_elapsed,
        time_per_step=time_per_step
    )
