fig, ax = plt.subplots(nrows=1, ncols=1)
legend_handler = PyPlot.matplotlib.legend_handler
collections = PyPlot.matplotlib.collections
lines = PyPlot.matplotlib.lines


TABLEAU_COLORS = PyPlot.matplotlib.colors.TABLEAU_COLORS


l1, = ax.plot([0, 1], [0, 1], label="Line 1")
l2, = ax.plot([0, 1], [0, 0.9], label="Line 2")
l3, = ax.plot([0, 1], [0, 0.8], label="Line 3")
l4, = ax.plot([0, 1], [0, 0.7], label="Line 4")


# le = ax.legend()
# le.get_legend_handler_map()
# l.HandlerLine2D(l1)
# handles, labels = ax.get_legend_handles_labels()

# ax.legend(handles, labels, handler_map=Dict(l1=> PyPlot.matplotlib.legend_handler.HandlerColorLineCollection(cmap="jet", numpoints=30)))
ax.legend([(l1, l2, l3, l4),l1], ["multi-color", "ef"] , handler_map=Dict((l1, l2, l3, l4) => legend_handler.HandlerTuple(ndivide=4, pad=0)))