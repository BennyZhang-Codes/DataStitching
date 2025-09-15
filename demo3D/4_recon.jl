_, _, _, kspha = get_kspace(hoseq; Δt=1);

fig = plt_kspha(kspha*2π, dt_adc; width=30, height=9, fontsize_label=6)
# fig.savefig("$(path)/out/kspha.png", dpi=300, transparent=true, bbox_inches="tight", pad_inches=0.0)


