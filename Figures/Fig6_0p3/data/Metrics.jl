using MAT
import Statistics: quantile
using PyPlot, PyCall
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
mcolors = matplotlib.colors

path = "$(@__DIR__)/Figures/Fig6/";
out_path     = "$(path)/out" ; if ispath(out_path) == false mkpath(out_path) end

# datas = [
#          "0p3_500_r30_pha_snr10_admm_20_TV_0.5.mat",
#          "0p3_500_r30_snr10_admm_20_TV_0.5.mat",
#          "0p3_500_r30_pha_snr10_cgnr_10_L2_0.0.mat",
#          "0p3_500_r30_pha_snr10_cgnr_20_L2_0.0.mat",
#          ]
# subfolders = [
#         "0p3_500_r30_pha_snr10_admm_20_TV_0p5",
#         "0p3_500_r30_snr10_admm_20_TV_0p5",
#         "0p3_500_r30_pha_snr10_cgnr_10_L2_0p0",
#         "0p3_500_r30_pha_snr10_cgnr_20_L2_0p0",
#          ]
datas = [
            "0p3_500_r30_144_ss1_snr10_cgnr_20_L2_1.0e-9.mat",
            "0p3_500_r30_144_ss3_snr10_cgnr_20_L2_1.0e-9.mat",
            "0p3_500_r30_196_ss1_snr10_cgnr_20_L2_1.0e-9.mat",
            "0p3_500_r30_196_ss3_snr10_cgnr_20_L2_1.0e-9.mat",
            "0p3_500_r30_256_ss1_snr10_cgnr_20_L2_1.0e-9.mat",
            "0p3_500_r30_256_ss3_snr10_cgnr_20_L2_1.0e-9.mat",
            ]
subfolders = [
           "0p3_500_r30_144_ss1_snr10_cgnr_20_L2_1e-9.mat",
           "0p3_500_r30_144_ss3_snr10_cgnr_20_L2_1e-9.mat",
           "0p3_500_r30_196_ss1_snr10_cgnr_20_L2_1e-9.mat",
           "0p3_500_r30_196_ss3_snr10_cgnr_20_L2_1e-9.mat",
           "0p3_500_r30_256_ss1_snr10_cgnr_20_L2_1e-9.mat",
           "0p3_500_r30_256_ss3_snr10_cgnr_20_L2_1e-9.mat",
            ]
for idx in [1, 2, 3, 4, 5, 6]
# idx = 2
subfolder = subfolders[idx];
out_path  = "$(path)/data/out/$(subfolder)" ; if ispath(out_path) == false mkpath(out_path) end

data = datas[idx]
data_file = "$(path)/data/$(data)"
@info "data file: $(data_file)"
data = matread(data_file);
imgs          = data["imgs"];
labels        = data["labels"];

x_ref         = data["x_ref"];
headmask      = data["headmask"];
B0map         = data["b0"];

# figure_width       = 5/2.53999863
# figure_height      = 5/2.53999863
fig_width          = 5
fig_height         = 5
vmaxp              = 99;
vminp              = 0;
color_facecolor    = "#ffffff";
dpi                = 900;

open("$(out_path)/Metrics.txt", "w") do file
    println(file, "$(datas[idx])\n")
# imgs   = imgStitched
# labels = labelStitched
nImg   = length(labels)
# vmin, vmax = quantile(abs.(imgs)[:], vminp/100), quantile(abs.(imgs)[:], vmaxp/100)
vmin, vmax = 0, 1
for idx = 1 : nImg
    @info labels[idx] 
    SSIM  = HO_SSIM(x_ref, abs.(imgs[:,:,idx])) 
    NRMSE = HO_NRMSE(x_ref, abs.(imgs[:,:,idx]))
    println(file, "$(labels[idx]) \n  SSIM: $(round(SSIM, digits=3)) \t NRMSE: $(round(NRMSE, digits=3))\n")

end
end
end