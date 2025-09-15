brainphantom = BrainPhantom(prefix="brain3D", x=0.5, y=0.5, z=0.5)
objbrain     = brainphantom
ss           = 4
start_end    = [130, 270]
db0_type    = :fat
db0_file    = :B0 
csm_type    = :real_32cha 
csm_nCoil   = 32
verbose     = false    

obj = brain_hophantom3D(brainphantom; 
                        ss=ss,
                        start_end=start_end,
                        db0_type=db0_type, db0_file=db0_file, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, verbose=verbose)

plt_phantom(obj, :ρ;)
plt_phantom(obj, :csm_mag;)
plt_phantom(obj, :Δw)

B0 = true    # turn on B0
T2 = true    # turn on T2

obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # set Δw to 0 if B0=false
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

