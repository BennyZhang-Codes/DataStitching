using KomaHighOrder
println("Testing...")


brainphantom = BrainPhantom("brain3D724")

obj = brain_hophantom2D(brainphantom)
plot_phantom_map(obj, :ρ)
plot_phantom_map_csm(obj, :mag; view_2d=true)