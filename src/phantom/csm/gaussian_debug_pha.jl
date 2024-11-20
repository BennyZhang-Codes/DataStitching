using PyPlot, LinearAlgebra

nX = 150
nY = 150
nCoil =100
nRow = 10
nCol = 10
relative_radius = 5
verbose = false

if nRow === nothing || nCol === nothing
    nRow, nCol = get_factors(nCoil)
else
    @assert nRow * nCol == nCoil "nCoil should be equal to nRow * nCol"
end

if verbose
    @info "Gaussian grid sensitivity" nX=nX nY=nY nRow=nRow nCol=nCol relative_radius=relative_radius
end

rX = nX/nRow/2
rY = nY/nCol/2

# Create the X and Y meshes
centerX = collect(range(rX, nX-rX, nRow));
centerY = collect(range(rY, nY-rY, nCol));

X = collect(range(0.5,nX-0.5,nX)) .* ones(1, nY);
Y = ones(nX) .* collect(range(0.5, nY-0.5,nY))';

sigma = 3
r_x = relative_radius*nX/nRow/2/sigma
r_y = relative_radius*nY/nCol/2/sigma

out = zeros(ComplexF64, (nX, nY, nCoil));
# Compute the Gaussian distribution for each BX BY box
for bx in 1:nRow
    for by in 1:nCol
        x = centerX[bx]
        y = centerY[by]

        mag = exp.(-(((X .- x)/r_x).^2 + ((Y .- y)/r_y).^2) ./ 2)

        pha = angle.((X .-  x) .+ im*(Y .- y))
        # pha = sqrt.((X.-x).^2 .+ (Y.-y).^2)./minimum([nX/nRow/2, nY/nCol/2]).^2 .* π
        # pha = (bx+by) % 2 == 0 ? pha : -pha 
        real = cos.(pha) .* mag
        imag = sin.(pha) .* mag
        out[:, :, (bx-1)*nCol+by] = (real .+ imag*im)
        # out[:, :, (bx-1)*nCol+by] = mag #.* exp.(im.*(bx*π+by*π))
    end
end
norm = sqrt.(sum(abs.(out) .^ 2, dims=3));
out = out./ norm;

# plt_image(abs.(sqrt.(sum(out.^2; dims=3))[:,:,1]))
# plt_image(sqrt.(sum(abs.(out).^2; dims=3)[:,:,1]))
# plt_image(angle.(out[:,:, 55]))
# plt_image(angle.(out[:,:, 56]))

# plt_image(angle.(out[:,:, 55]).-angle.(out[:,:,56]), vmin=-π, vmax=π)
# plt_image(angle.(out[:,:, 55]).-angle.(out[:,:,65]), vmin=-π, vmax=π)

fig = plt_images(abs.(out); dim=3, nRow=nRow, nCol=nCol)
fig = plt_images(angle.(out); dim=3, nRow=nRow, nCol=nCol)
fig.savefig("E:/skope/20241213/csm_gaussian_grid_block/gaussian_pha.png", dpi=300, transparent=true, bbox_inches="tight", pad_inches=0)

# fig, ax = plt.subplots(1, 1)
# for row = 1:nRow
#     ax.plot(abs.(out[7, :, row]))
# end

fig, axs = plt.subplots(1, 2)
ax1 = axs[1]
ax2 = axs[2]
for col = 1:nCol
    ax1.plot(abs.(out[nX÷2, :, (nRow÷2-1)*nCol+col]))
end
ax2.scatter(collect(1:1:nCoil), abs.(out)[75, 75, :], s=5)
ax2.scatter(collect(1:1:nCoil), angle.(out)[75, 75, :], s=5)



# fig, axs = plt.subplots(1, 2)
# ax1 = axs[1]
# ax2 = axs[2]
# for row = 1:nRow
#     for col = 1:nCol
#         if row == col
#             ax1.plot(abs.(diag(out[:, :, (row-1)*nCol+col])))
#         end
#     end
# end
# ax2.set_ylim(-π, π)
# for row = 1:nRow:nX
#     for col = 1:nCol:nY
#         ax2.scatter(ones(nCoil)*(row*nY+col), angle.(out)[row, col, :], s=10)
#     end
# end
