import KomaMRI.KomaMRIPlots: interp_map, theme_chooser, plot_koma
import KomaMRI.KomaMRIPlots: plot_kspace, plot_seqd, plot_seq, plot_phantom_map

include("plotly/plotly.jl")

include("plt/plt.jl")

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
		plot_bgcolor = "rgb(242,242,242)"
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
