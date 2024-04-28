include("plot_kspace.jl")
include("plot_seq.jl")
include("plot_seqd.jl")
include("plot_traj2d.jl")

export plot_seq
export plot_seqd, plot_hoseqd
export plot_grads_cumtrapz
export plot_traj2d

function HO_theme_chooser(mode::Symbol)
    @assert mode in [:light, :light_fontblack, :dark, :white] "the mode must be :light, :light_fontblack, :dark, or :white"
	if mode == :dark
		bgcolor = "rgba(0,0,0,0)"#"rgb(13,16,17)"
		text_color = "gray"
		plot_bgcolor = "rgb(22,26,29)" #"rgb(33,37,41)"
		grid_color = "rgb(40,52,66)" #rgb(40,40,40)
		sep_color = "white"
    elseif mode == :light
		bgcolor = "rgba(0,0,0,0)"#"white"
		text_color = "gray"#"rgb(49,70,101)"
		plot_bgcolor = "rgb(229,236,246)"
		grid_color = "white"
		sep_color = "black"
    elseif mode == :light_fontblack
        bgcolor = "rgba(0,0,0,0)"#"white"
		text_color = "black"#"rgb(49,70,101)"
		plot_bgcolor = "rgb(229,236,246)"
		grid_color = "white"
		sep_color = "black"
	elseif mode == :white
        bgcolor = "rgba(0,0,0,0)"#"white"
		text_color = "black"#"rgb(49,70,101)"
		plot_bgcolor = "white"
		grid_color = "rgb(40,52,66)"
		sep_color = "black"
    end
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color
end



"""
hoseq = demo_hoseq()
thememode = :white

plot_grads_cumtrapz(hoseq; thememode=thememode)
plot_grads_cumtrapz(hoseq, 2; thememode=thememode)
plot_kspace(hoseq; thememode=thememode)
plot_kspace(hoseq, :x; thememode=thememode)

plot_seq(hoseq; thememode=thememode)

plot_seqd(hoseq; thememode=thememode)
plot_hoseqd(hoseq; thememode=thememode)
"""

