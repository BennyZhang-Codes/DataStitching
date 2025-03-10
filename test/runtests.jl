using HighOrderMRI
println("Testing...")


brainphantom = BrainPhantom("brain3D724")


# setting the coil sensitivity used in the simulation
csm_type  = :birdcage;      # a simulated birdcage coil-sensitivity
csm_nCoil = 9;              # 8-channel
csm_nRow  = 3;
csm_nCol  = 3;

obj = brain_hophantom2D(brainphantom; csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol)
plot_phantom_map(obj, :ρ; view_2d=true)
plot_phantom_map_csm(obj, :mag; coil_idx=1, view_2d=true)
plot_phantom_map_csm(obj, :pha; coil_idx=1, view_2d=true)

# obj = brain_phantom2D_reference(brainphantom; key=:csm, target_fov=(150,150), target_resolution=(1,1),
#     csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol)