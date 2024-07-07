b2d = BrainPhantom()
b3d_171 = brain3D(file = phantom_dict[:brain3d_171])
b3d_285 = brain3D(file = phantom_dict[:brain3d_285])

obj2d = brain_phantom2D(b2d); info(obj2d);
obj3d_171 = brain_phantom3D(b3d_171;ss=3, start_end=[180,220]); info(obj3d_171);
obj3d_285 = brain_phantom3D(b3d_285;ss=3, start_end=[180,220]); info(obj3d_285);

# obj3d_171 = brain_phantom3D(b3d_171;ss=1, start_end=[1,400]); info(obj3d_171);
# obj3d_285 = brain_phantom3D(b3d_285;ss=1, start_end=[1,400]); info(obj3d_285);

obj.Δw .= obj.Δw * 0;

plot_phantom_map(obj2d, :ρ);
plot_phantom_map(obj3d_171, :ρ);
plot_phantom_map(obj3d_285, :ρ);


