using KomaHighOrder
import KomaMRI.KomaMRIBase: get_samples

import KomaMRI.KomaMRIPlots: plot_seq

import KomaMRI.KomaMRIPlots: plot_seqd

import KomaMRI.KomaMRICore: run_spin_precession!, run_spin_excitation!
import KomaMRI.KomaMRICore: output_Ndim, initialize_spins_state, Bloch, sim_output_dim

import KomaMRI.KomaMRICore: run_spin_excitation_parallel!, run_spin_precession_parallel!, run_sim_time_iter!

import KomaMRI.KomaMRIBase: discretize, get_theo_Gi, get_grads

import KomaMRI.KomaMRICore: simulate





methods(KomaHighOrder.simulate)
methods(KomaHighOrder.run_sim_time_iter!)
methods(KomaHighOrder.run_spin_precession_parallel!)
methods(KomaHighOrder.run_spin_excitation_parallel!)

methods(KomaHighOrder.run_spin_precession!)
methods(KomaHighOrder.run_spin_excitation!)
methods(KomaHighOrder.initialize_spins_state)
methods(KomaHighOrder.output_Ndim)
methods(KomaHighOrder.sim_output_dim)

methods(KomaHighOrder.discretize)
methods(KomaHighOrder.get_theo_Gi)
methods(KomaHighOrder.get_grads)
methods(KomaHighOrder.get_samples)

