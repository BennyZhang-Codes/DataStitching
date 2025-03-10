import KomaMRI.KomaMRICore: run_spin_precession!, run_spin_excitation!
import KomaMRI.KomaMRICore: output_Ndim, initialize_spins_state, Bloch, sim_output_dim

import KomaMRI.KomaMRICore: mul!, Q
import KomaMRI.KomaMRICore: run_spin_excitation_parallel!, run_spin_precession_parallel!, run_sim_time_iter!
import KomaMRI.KomaMRICore: run_spin_precession!, run_spin_excitation!
import KomaMRI.KomaMRICore: simulate

# define BlochHighOrder struct first, for workflow control in the run_spin_excitation! & run_spin_precession! functions
include("BlochHighOrder.jl")
export BlochHighOrder

# define simulation functions for different phantoms
include("Phantom/simulation.jl")       # simulation for Phantom defined in KomaMRI.jl
include("HO_Phantom/simulation.jl")    # simulation for HO_Phantom defined in HighOrderMRI.jl

include("SimType.jl")
export SimType