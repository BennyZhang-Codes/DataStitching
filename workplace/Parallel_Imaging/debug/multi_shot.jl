using KomaHighOrder, PyPlot

outpath = "$(@__DIR__)/workplace/Parallel_Imaging/debug/out/multi_shot_0.1"; if ispath(outpath) == false mkpath(outpath) end     # output directory
fileprefix = "spiral_R30"
##############################################################################################
# Setup
##############################################################################################
T = Float32;
R = 30
matrix_origin = 500
matrix_target = 500
Nx = Ny = matrix_target
grad_scale    = 1/(matrix_origin/matrix_target)



seq = load_seq(seqname="spiral", r=R)[4:end]
hoseq = HO_Sequence(seq)

seq_delayTR = Sequence([seq.GR[3, end] Grad(0,0); seq.GR[3, end] Grad(0,0); Grad(0,10) Grad(0,0)])


gx = seq.GR[1, end-1] * grad_scale
gy = seq.GR[2, end-1] * grad_scale

nInterleaves = 5
⦞ = collect(range(0,2,nInterleaves+1))[1:end-1]

seqs = []
seq_ms = Sequence()
for ∠ in ⦞
    g1 = rotz(∠*π) * transpose([gx.A gy.A gx.A*0])
    seq1 = copy(seq)
    seq1.GR[1, end-1] = Grad(g1[1, :], gx.T, gx.rise, gx.fall, gx.delay)
    seq1.GR[2, end-1] = Grad(g1[2, :], gy.T, gy.rise, gy.fall, gy.delay)
    seq_ms += seq1 + seq_delayTR
    push!(seqs, seq1)
end

# plot_seq(seq_ms)
# plot_kspace(seq_ms)
##### 1. sequence
hoseq = HO_Sequence(seq_ms)

_, k_nominal, _, _ = get_kspace(hoseq; Δt=1);
fig_traj = plt_traj(k_nominal'; color_label="#CCCCCC", width=8, height=8)
fig_traj.savefig("$(path)/$(fileprefix)_$(nInterleaves)shot_traj.png", dpi=300, bbox_inches="tight", pad_inches=0, transparent=true)



# sys = Scanner();
# sim_params = KomaMRICore.default_sim_params(); 
# sim_params["sim_method"] = BHO;
# sim_params["gpu"] = true;
# sim_params["return_type"]="state";
# sim_params["precision"] = "f64"
# sim_params["Nblocks"] = 1000
# # # 4. simulate
# state = simulate(obj, HO_Sequence(seq_ms), sys; sim_params);
# p_xy_mag = plot_mag(obj, state, :xy_mag; title="Mxy_mag", darkmode=true, view_2d=true)
# p_xy_pha = plot_mag(obj, state, :xy_pha; title="Mxy_pha", darkmode=true, view_2d=true)
# p_z      = plot_mag(obj, state, :z     ; title="Mz"     , darkmode=true, view_2d=true)





simtype  = SimType(B0=false, T2=false, ss=1)
BHO = BlochHighOrder("000", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.1, y=0.1, z=0.2) # decide which phantom file to use
##### 2. phantom
# setting the coil sensitivity used in the simulation
csm_type  = :fan;      # a simulated birdcage coil-sensitivity
csm_nCoil = 1;         # 1-channel
csm_nRow  = 1;
csm_nCol  = 1;

db0_type  = :quadratic;     
db0_max   = :0.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map

obj = brain_hophantom2D(phantom; ss=simtype.ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
                        db0_type=db0_type, db0_max=db0_max); obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; # TODO: fix the bug: gre 

##### 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BHO;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";
sim_params["precision"]   = "f64"
sim_params["Nblocks"] = 1000


signal = Array{Complex{T}}(undef, 21000*length(seqs), 1, 1);
for idx in eachindex(seqs)
    s = seqs[idx]
    signal[1 + (idx-1)*21000:idx*21000, :, :] = simulate(obj, HO_Sequence(s), sys; sim_params);
end

raw = signal_to_raw_data(signal, HO_Sequence(seq_ms), :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);

fig = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig.savefig("$(path)/$(fileprefix)_$(nInterleaves)shot_nufft.png", dpi=300, bbox_inches="tight", pad_inches=0, transparent=true)


##### 4. simulate
# signal = simulate(obj, HO_Sequence(seq_ms), sys; sim_params); figure(); plt.plot(signal[:])
# raw = signal_to_raw_data(signal, HO_Sequence(seq_ms), :nominal; sim_params=copy(sim_params));
# img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
# fig = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]); width=12/2.54, height=12/2.54, title="Sim: $(BHO.name)")

# fig = plt_images(permutedims(mapslices(rotl90, img_nufft,dims=[1,2]), [3, 1, 2]),width=10, height=5)

matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=100)
# matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 8/2.54
figure_height      = 4/2.54
linewidth          = 0.5
ticklength         = 1.5
fontsize_legend    = 5
fontsize_label     = 6
fontsize_ticklabel = 4
fontsize_subfigure = 8
pad_labeltick      = 2
pad_label          = 2
color_facecolor    = "#1f1f1f"
color_label        = "#cccccc"
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor)

ax.tick_params(axis="both", length=ticklength, width=linewidth, color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_color(color_label)
    ax.spines[spine].set_visible(false)
end
ax.set_facecolor(color_facecolor)
ax.plot(real.(signal[:]))
ax.set_ylabel("Signal [a.u.]", fontsize=fontsize_label, color=color_label)
# ax.set_xlabel("Samples", fontsize=fontsize_label, color=color_label)
fig.tight_layout()
fig.savefig("$(outpath)/$(fileprefix)_$(nInterleaves)shot_signal.png", dpi=300, bbox_inches="tight", pad_inches=0, transparent=true)

