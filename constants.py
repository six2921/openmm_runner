from openmm import unit
from openmm.app import HBonds

# Integrator settings
integrator_config = {
    "temp": 300 * unit.kelvin,
    "dt": 0.002 * unit.picoseconds,
    "gamma": 1.0 / unit.picosecond
}

# System settings
system_config = {
    "pressure": 1.0 * unit.atmospheres,
    "baro_temp": 300 * unit.kelvin,
    "baro_freq": 25
}

# Forcefield settings
forcefield_config = {
    "forcefields": ["amber/ff14SB.xml", "amber/tip3p_standard.xml"],
    "forcefield_kwargs": {
        "constraints": HBonds,
        "rigidWater": True,
        "removeCMMotion": False,
    },
    "small_molecule_ff": "gaff-2.11"
}
