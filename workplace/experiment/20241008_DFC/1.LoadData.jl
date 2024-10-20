

path = "/home/jyzhang/Desktop/pulseq/20241010_skope_fa90/invivo/"
seq_file = "$(path)/seq/xw_sp2d_7T-1mm-200-r4-noSync-fa90.seq"
dfc_file = "$(path)/dfc/xw_sp2d_7T-1mm-200-r4.mat"

# 1. read the *.seq file and reverse the sign of the gradient (axis x)
seq = read_seq(seq_file)[end-9:end-3]; 
seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)

# 2. load the *.mat file and extract the DFC data, then, composing the HO_Sequence object
dfcStitched, dfcStandard, ntStitched, ntStandard = load_dfc_mat(dfc_file);
hoseqStitched = HO_Sequence(seq); # hoseq, defined in KomaHighOrder.jl
hoseqStandard = HO_Sequence(seq); # hoseq, defined in KomaHighOrder.jl
hoseqStitched.GR_dfc[2:4, :] = hoseqStitched.SEQ.GR; # copy the 1st-order gradient data from the seq object to the hoseq object
hoseqStandard.GR_dfc[2:4, :] = hoseqStandard.SEQ.GR; # copy the 1st-order gradient data from the seq object to the hoseq object
hoseqStitched.GR_dfc[:,6] = dfcStitched;    # "6" is the index of the readout block in the spiral sequence
hoseqStandard.GR_dfc[:,6] = dfcStandard;    # "6" is the index of the readout block in the spiral sequence
plot_seq(hoseqStitched)
plot_seq(hoseqStandard)
# finally, hoseq* contains both the nominal trajectory and the measured trajectory (up to 2nd-order)


_, K_nominal, _, K_dfcStitched = get_kspace(hoseqStitched; Δt=1);
_, _, _, K_dfcStandard = get_kspace(hoseqStandard; Δt=1);

times = KomaMRIBase.get_adc_sampling_times(seq);

tr_nominal          = Trajectory( K_nominal'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfcStitched     = Trajectory(K_dfcStitched'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfcStandard     = Trajectory(K_dfcStandard'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);


