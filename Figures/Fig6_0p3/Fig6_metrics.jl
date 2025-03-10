include("$(@__DIR__)/Figures/plt_preset.jl");
using MAT
import Statistics: quantile

outpath = "$(@__DIR__)/Figures/Fig6/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory

###### load images reconstructed by different considerations of the extended signal model 
solver = "admm"; regularization = "TV"; λ = 1.e-4; iter=20;
solver = "cgnr"; regularization = "L2"; λ = 1.e-3; iter=20;
snr=10;
matfile  = "R30_snr$(snr)_$(solver)_$(iter)_$(regularization)_$(λ)"

matdatas = MAT.matread("$(outpath)/$(matfile).mat")
imgs     = matdatas["imgs"];
labels   = matdatas["labels"];
B0map    = matdatas["B0map"];
headmask = matdatas["headmask"];
x_ref    = matdatas["x_ref"];
nFrame, nX, nY = size(imgs)


###### calculate metrics
for i in 1:nFrame
    img = imgs[i, :, :];
    # plt_image(img)
    # label = labels[i];
    ssim = HO_SSIM(img,x_ref)
    # println(label, " SSIM: ", ssim)
    nrmse = HO_NRMSE(x_ref, img)
    # println(label, " NRMSE: ", nrmse)
    println(ssim)
end
img = imgs[3, :, :]
HO_SSIM(x_ref,x_ref)

vmaxp              = 99;
vminp              = 1;
vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
fig = plt_image(abs.(imgs[3, :, :] - imgs[4, :, :]) * 5, vmin=vmin, vmax=vmax)
fig.savefig("$(outpath)/Fig6_GT_$(matfile)_diff.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0)


plt_image(abs.(imgs[3, :, :] - HO_img_scale(imgs[3, :, :], x_ref)))
plt_image(abs.(imgs[4, :, :] - HO_img_scale(imgs[4, :, :], x_ref)))