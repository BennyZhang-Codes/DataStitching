using HighOrderMRI

brainphantom = BrainPhantom("brain3D724")
# setting the coil sensitivity used in the simulation
csm_type  = :rect_gaussian;      # a simulated birdcage coil-sensitivity
csm_nCoil = 9;              # 8-channel
csm_nRow  = 3;
csm_nCol  = 3;

db0_type  = :quadratic;     
db0_max   = :125.;

obj = brain_hophantom2D(brainphantom; 
                        db0_type=db0_type, db0_max=db0_max, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol)
plot_phantom_map(obj, :ρ; view_2d=true)
plot_phantom_map_csm(obj, :mag; coil_idx=5, view_2d=true)
plot_phantom_map_csm(obj, :pha; coil_idx=5, view_2d=true)


brain_phantom2D_reference(brainphantom, :csm, (150.,150.), (1.,1.);
                        db0_type=db0_type, db0_max=db0_max, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol)