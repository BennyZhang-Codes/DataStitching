using MAT, NIfTI


path = "/home/jyzhang/Desktop/pulseq/20241104_ABDL/contour"
mat_file = "$(path)/gre_meas_MID00117_FID53005_pulseq_v0_gres6_1p0_standard.mat"
imgs = matread(mat_file)["gre_imgs"]
imgs = mapslices(rotr90, imgs, dims=[2,3])
raw  = permutedims(imgs, (2, 3, 1))[:, :, :]

ni = niread("T2-flair.nii")
# ni.header.dim = (3, 200, 200, 6, 1, 1, 1, 1)


# nii = NIVolume(ni.raw; voxel_size=(1,1,2), dim_info=(1,2,3))
nii = NIVolume(ni.raw)
niwrite("$(path)/gres6_1p0.nii", nii)


nCha, nX, nY = size(imgs)
for cha = 1:nCha
    nii = NIVolume(permutedims(imgs, (2, 3, 1))[:, :, cha:cha])
    niwrite("$(path)/gres6_1p0_echo$(cha).nii", nii)
end


#= Bash command
mri_synthseg --i gres6_1p0.nii.gz --o gres6_1p0_synthseg.nii.gz --threads 40 --robust --resample gres6_1p0_resample.nii.gz
hd-bet -i T2-flair.nii -o T2-flair_brain.nii.gz
=#