using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISimulation
using PlotlyJS, MAT, ImageQualityIndexes, ImageDistances
using ProgressMeter

############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=true, T2=false, ss=5)
BHO = BlochHighOrder("000", true, false)
dir = "$(@__DIR__)/B0/results/spinstate"; if ispath(dir) == false mkpath(dir) end
maxOffresonance = 0.
hoseq_base = demo_hoseq();# plot_hoseqd(hoseq);
key = :fatsat
hoseq_dict = Dict(
    :fatsat      =>hoseq_base[2:3],
    :fatsat_excit=>hoseq_base[3:5],
    :excit       =>hoseq_base[4:5],
    )

# for key in keys(hoseq_dict)
    prefix = "dB0precessionOff_dB0max$(maxOffresonance)"
    prefix = BHO.Δw_excitation ? "$(prefix)_dB0excitOn_$(String(key))" : "$(prefix)_dB0excitOff_$(String(key))"

    ############################################################################################## 
    # Simu
    ############################################################################################## 
    # 1. hoseq

    hoseq = hoseq_dict[key]
    p_hoseq = plot_seq(hoseq)

    # 2. phantom
    obj = brain_hophantom2D(BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2); ss=simtype.ss, location=0.8, B0type=:quadratic, maxOffresonance=maxOffresonance)
    obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
    obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation
    # 3. scanner & sim_params
    sys = Scanner();
    sim_params = KomaMRICore.default_sim_params()
    sim_params["sim_method"] = BHO;
    sim_params["gpu"] = true;
    sim_params["return_type"]="state";
    sim_params["precision"] = "f64"
    sim_params["Nblocks"] = 1
    # 4. simulate
    state = simulate(obj, hoseq, sys; sim_params);
    p_xy_mag = plot_mag(obj, state, :xy_mag; darkmode=true, view_2d=true)
    p_xy_pha = plot_mag(obj, state, :xy_pha; darkmode=true, view_2d=true)
    p_z      = plot_mag(obj, state, :z     ; darkmode=true, view_2d=true)

    savefig(p_hoseq , "$(dir)/$(prefix)_hoseq.svg" , width=500,height=500,format="svg")
    savefig(p_xy_mag, "$(dir)/$(prefix)_xy_mag.svg", width=500,height=500,format="svg")
    savefig(p_xy_pha, "$(dir)/$(prefix)_xy_pha.svg", width=500,height=500,format="svg")
    savefig(p_z     , "$(dir)/$(prefix)_z.svg"     , width=500,height=500,format="svg")
# end
