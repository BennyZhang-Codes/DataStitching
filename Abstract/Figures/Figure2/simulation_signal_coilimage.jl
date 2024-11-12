using KomaHighOrder
using MAT, PyPlot
import Statistics: quantile
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=true, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
csmtype= :birdcage; nCoil=8; nrows=2; ncols=4;
BHO = BlochHighOrder("000", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
maxOffresonance = 100.                                           # set maximum off-resonance frequency in Hz for quadratic B0 map
Nx = Ny = 150;

T = Float64;

outpath = "workplace/Abstract/Figures/Figure2/out"; if ispath(outpath) == false mkpath(outpath) end


#############################################################################
# 1. Simulation
#############################################################################
# 1. sequence
hoseq_stitched = load_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseq_standard = load_hoseq(dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, csmtype=csmtype, nCoil=nCoil, B0type=:quadratic, maxOffresonance=maxOffresonance); 
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"
sim_params["gpu_device"] = 1
# 4. simulate
signal = simulate(obj, hoseq_stitched, sys; sim_params);
# data = signal[:,:,1];
# raw = signal_to_raw_data(signal, hoseq_stitched, :nominal; sim_params=copy(sim_params));
# img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
# fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]); width=12/2.54, height=12/2.54)
# fig_cha = plt_images(permutedims(mapslices(rotl90, img_nufft,dims=[1,2]), [3, 1, 2]),width=8, height=8)


#############################################################################
# 2. Adding noise to signal data
#############################################################################
MAT.matwrite(outpath *"/fully_snr10_db0100.mat", Dict("signal"=>signal))



signal = matread("$(outpath)/fully_snr10_db0100.mat")["signal"];

snr = 10;

data = signal[:,:,1];
nSample, nCha = size(data);
signalAmpl = sum(abs.(data), dims=1)/ nSample;
noise = signalAmpl/snr .* ( randn(size(data))+ 1im*randn(size(data)));
data = data + noise;
# plt.plot(abs.(data[:, 1]), linewidth=0.5)
raw = signal_to_raw_data(reshape(data, (nSample, nCha, 1)), hoseq_stitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]); width=12/2.54, height=12/2.54)

imgs_coil = permutedims(mapslices(rotl90, img_nufft,dims=[1,2]), [3, 1, 2]);
fig_cha = plt_images(imgs_coil,width=8, height=8)




nSample, nCha = size(data);
###############################
# Plotting signal for each channel
###############################
matplotlib.rc("mathtext", default="regular")
# matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rcParams["mathtext.default"]
figure_width       = 3.5/2.54
figure_height      = 1.8/2.54
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

for cha = 1 : nCha
    n = real.(noise[:, cha]);
    s = abs.(data[:, cha]);
    ymax = maximum(s);

    # fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor)
    # ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    # for spine in ax.spines  # "left", "right", "bottom", "top"
    #     ax.spines[spine].set_color(color_label)
    #     ax.spines[spine].set_visible(false)
    # end
    # ax.set_facecolor(color_facecolor)
    # ax.set_ylim(0, ymax)
    # ax.plot(s, linewidth=0.5, color="C$(cha%9-1)")
    # fig.tight_layout(pad=0)
    # fig.savefig("$(outpath)/Fig_signal_cha$(cha).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
    # fig.savefig("$(outpath)/Fig_signal_cha$(cha).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)

    figure_width       = 2/2.54
    figure_height      = 1/2.54
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecolor)
    ax.set_ylim(minimum(real.(noise[:, cha])), ymax/20)
    ax.plot(n, linewidth=0.05, color="C$(cha%9-1)")
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/Fig_noise_cha$(cha).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
    fig.savefig("$(outpath)/Fig_noise_cha$(cha).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
end

###############################
# Plotting signal for each channel
###############################
figure_width       = 5/2.54
figure_height      = 5/2.54

imgs = imgs_coil;
vmin, vmax = quantile(abs.(imgs)[:], 0), quantile(abs.(imgs)[:], 0.99999)
for idx = 1 : nCha
    img = abs.(imgs)[idx, :,:]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/Fig_coilimage_cha$(idx).png", dpi=900, transparent=true, bbox_inches="tight", pad_inches=0)
end




