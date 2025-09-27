using Test
using HighOrderMRI

brainphantom = BrainPhantom(prefix="brain3D", x=0.5, y=0.5, z=0.5)

objbrain     = brainphantom
axis_SAG     = "axial"     
ss           = 4
location     = [100, 120]        
b0_type      = :real
csm_type     = :real_32cha 
csm_nCoil    = 32
verbose      = false    

obj = brain_hophantom3D(brainphantom; 
                        axis=axis_SAG, 
                        ss=ss, 
                        location=location, 
                        b0_type=b0_type, 
                        csm_type=csm_type, 
                        csm_nCoil=csm_nCoil, 
                        verbose=verbose)

# plt_phantom(obj, :ρ; view_2d=false)

nSpin, nCoil = size(obj.csm)
@test nSpin == length(obj.x)
@test nCoil == csm_nCoil


