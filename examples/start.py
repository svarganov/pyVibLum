import numpy as np
from pathlib import Path
import pyviblum.viblum as vl


# Set path to input data. 
base_path = Path.cwd()
subdir = "data"

path_to_elec  = base_path.joinpath(subdir,"OpenMolcas.out")
path_to_grads = base_path.joinpath(subdir,"gradients")
path_to_nacs  = base_path.joinpath(subdir,"nacs")
path_to_coord = base_path.joinpath(subdir,"GCoord.xyz")
path_to_hess  = base_path.joinpath(subdir,"GHessian.dat")

path_to_data = {
   "electronic":  path_to_elec,
   "gradients":   path_to_grads,
   "couplings":   path_to_nacs,
   "coordinates": path_to_coord, 
   "hessian":     path_to_hess
}


# Electronic (CASSCF) states included in the spin-orbit coupling calculations.
states_electronic = list(range(1, 22, 1))


# Spin-orbit states to be included in luminescence calculations. The last state is the emitting state.
states_spin_orbit = list(range(1, 17, 2))
states_spin_orbit.append(53)


# Spin multiplicity.
S = 1.5


# Hessian format.
hessian_format = "gamess"


# Simulation parameters.
sim_params = {
   "line_width":  3.5,
   "min_energy":  2.5950e0,
   "max_energy":  2.7164e0,
   "step":        0.0001e0,
   "num_exct":    3,
   "num_modes":   6
}


# Wrap input data and simulation parameters.
input_and_sim_params = {
    "paths": path_to_data,
    "elstates": states_electronic,
    "sostates": states_spin_orbit,
    "spin": S,
    "format": hessian_format,
    "params": sim_params
}


# Calculate vibronic intensities.
vl.collect_data_for_simulation(**input_and_sim_params)

