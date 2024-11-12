using PyPlot, MAT
using KomaHighOrder
sh = SphericalHarmonics()

path         = "E:/pulseq/20241010_skope_fa90/invivo"
r=4
seq_file     = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin("r$(r)", f)][1]
dfc_file     = "$(path)/dfc/" * [f for f in readdir("$(path)/dfc") if occursin("r$(r).mat", f)][1]

outpath = "workplace/Abstract/Figures/Figure4/out"; if ispath(outpath) == false mkpath(outpath) end


########################################################################
# 1. load the *.seq file and extract some infomations
########################################################################
seq = read_seq(seq_file)[end-9:end-3]; 
seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
t_gradFreeDelay = seq.DEF["skope_gradFreeDelay"]; # [s]


########################################################################
# 2. load the *.mat file and extract the DFC data
########################################################################

dt            = matread(dfc_file)["dt"];  # [s]
dfc           = matread(dfc_file)["bfieldStitched"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
nSampleAllSeg = matread(dfc_file)["nSampleAllSegStitched"]; 
nGradPoint, nTerm = size(matread(dfc_file)["bfieldStitched"]);   

t = dt * (nGradPoint-1);
GR_dfc = reshape([KomaMRIBase.Grad(dfc[idx,:], t, dt/2, dt/2, 0) for idx=1:9], :, 1);

hoseq = HO_Sequence(seq);                     # hoseq, defined in KomaHighOrder.jl
hoseq.GR_dfc[2:4, :] = hoseq.SEQ.GR;  # copy the 1st-order gradient data from the seq object to the hoseq object
hoseq.GR_dfc[:,6] = GR_dfc;              # "6" is the index of the readout block in the spiral sequence
plot_seq(hoseq)
# finally, hoseq* contains both the nominal trajectory and the measured trajectory (up to 2nd-order)




samples = get_samples(hoseq; off_val=Inf);
t_adc = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);


a = Int64.(vec(nSampleAllSeg)); a[1] += 3; e = cumsum(a); s = e .- a;
t_trigger = (s .+ 1)*dt .+ KomaMRIBase.get_block_start_times(hoseq.SEQ)[5] .+ t_gradFreeDelay;


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 9/2.54
figure_height      = 3/2.54
linewidth          = 0.5
ticklength         = 1.5
fontsize_legend    = 5
fontsize_label     = 6
fontsize_ticklabel = 4
fontsize_subfigure = 8
pad_labeltick      = 2
pad_label          = 2
color_facecolor    = "#ffffff"
color_label        = "#000000"
color_trig         = "#234765"
color_trig         = "#000000"

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)

ax = axs[1, 1]
ax.tick_params(axis="both", length=ticklength, width=linewidth, color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_color(color_label)
    ax.spines[spine].set_visible(false)
end
ax.set_xlim(samples.gz.t[1]*1e3, samples.gz.t[end]*1e3)
ax.set_facecolor(color_facecolor)
ax.yaxis.set_major_locator(plt.MultipleLocator(20))

ax.axhline(y=0, color="#555555", linewidth=linewidth)

ax.plot(samples.gx.t*1e3,       samples.gx.A*1e3, color="#636EFA", linewidth=linewidth, label=L"G_{x}")
ax.plot(samples.gy.t*1e3,       samples.gy.A*1e3, color="#EF553B", linewidth=linewidth, label=L"G_{y}")
ax.plot(samples.gz.t*1e3,       samples.gz.A*1e3, color="#00CC96", linewidth=linewidth, label=L"G_{z}")

ax.plot(samples.rf.t*1e3, abs.(samples.rf.A)*1e6, color="#AB63FA", linewidth=linewidth, label=L"|B_{1}|")

ax.scatter(t_trigger*1e3, -55*ones(length(t_trigger)), s=15, marker=L"\mapsup", color=color_trig, linewidth=0.1, label=L"Stitching \ trigger")
ax.scatter(t_trigger[1]*1e3, -40, s=15, marker=L"\twoheaduparrow", color=color_trig, linewidth=0.1, label=L"Standard \ trigger")

ax.set_ylim(-70, 70)
ax.set_ylabel("Amplitude (mT/m, μT)", fontsize=fontsize_label, color=color_label)
ax.set_xlabel("Time (ms)", fontsize=fontsize_label, color=color_label)
ax.legend(loc="upper left", bbox_to_anchor=(0, 1.1), fontsize=fontsize_legend, labelcolor=color_label, ncols=6, frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1)

fig.tight_layout(pad=0, h_pad=0.5, w_pad=0.5)
fig.savefig("$(outpath)/Fig_sequence_r$(r).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
fig.savefig("$(outpath)/Fig_sequence_r$(r).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)