brainphantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2)
# setting the coil sensitivity used in the simulation
csm_type  = :rect_gaussian;      # a simulated birdcage coil-sensitivity
csm_nCoil = 9;              # 8-channel
csm_nRow  = 3;
csm_nCol  = 3;

db0_type  = :quadratic;     
db0_max   = :125.;

obj = brain_hophantom2D(brainphantom; 
                        ss=10,
                        db0_type=db0_type, db0_max=db0_max, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol)
# plt_phantom_map(obj, :ρ; view_2d=true)

nSpin, nCoil = size(obj.csm)

@test nSpin == length(obj.x)
@test nCoil == csm_nCoil
