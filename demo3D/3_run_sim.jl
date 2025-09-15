BHO = BlochHighOrder("1111", true, true)            

# 2. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params();
sim_params["sim_method"]  = BHO;      # using "BlochHighOrder" for simulation with high-order terms
sim_params["return_type"] = "mat";    # setting with "mat", return the signal data for all channel
sim_params["precision"]   = "f64";
sim_params["gpu"]         = true;     # using GPU for simulation
sim_params["gpu_device"]  = 0;        # set the GPU device number, if using GPU for simulation
sim_params["Nblocks"]     = 2000;     # the number of blocks for GPU simulation, set according to the GPU memory

# 3. simulate
signal    = simulate(obj, hoseq, sys; sim_params);          


s = reshape(signal, nSample_adc, :, csm_nCoil); # reshape to [nSample_adc, nADC, nCha]
plt_plot([abs.(s[:,1,8])])


fig = plt_plot([abs.(signal[:,1])]; width=30, height=5, fontsize_label=6, ylabel="Signal Amplitude [a.u.]", xlabel="ADC Sample")
# fig.savefig("$(path)/out/signal.png", dpi=300, transparent=true, bbox_inches="tight", pad_inches=0.0)



