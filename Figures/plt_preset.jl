using PyPlot, PyCall
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

mtransforms = PyPlot.plt.matplotlib.transforms
mcolors = matplotlib.colors

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["axes.titlepad"] = 2

matplotlib.rc("mathtext", default="regular");
matplotlib.rc("figure", dpi=200);
matplotlib.rc("font", family="Times New Roman");
matplotlib.rcParams["mathtext.default"];

figure_width       = 5/2.54;
figure_height      = 4/2.54;

vmaxp              = 99;
vminp              = 1;
cmap               = "gray";

fontsize_legend    = 5;
fontsize_label     = 6;
fontsize_subfigure = 8;
fontsize_ticklabel = 4;

linewidth          = 0.5;
ticklength         = 1.5;

pad_label          = 2;
pad_labeltick      = 2;

color_facecolor    = "#ffffff";
color_label        = "#000000";
color_subfigure    = "#ffffff";



