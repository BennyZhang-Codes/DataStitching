using KomaHighOrder
using PyPlot, LinearAlgebra

nX    = 150
nY    = 150
nRow = 5
nCol = 5
nCoil = nRow*nCol
nBlock = 4

nRowall  = nRow * nBlock
nColall  = nCol * nBlock

relative_radius = 5
verbose = false

if nRow === nothing || nCol === nothing
    nRow, nCol = get_factors(nCoil)
else
    @assert nRow * nCol == nCoil "nCoil should be equal to nRow * nCol"
end

if verbose
    @info "Gaussian grid sensitivity" nX=nX nY=nY nRow=nRow nCol=nCol nBlock=nBlock nCoil=nCoil nCoil_all=nRowall*nColall relative_radius=relative_radius 
end

rX = nX/nRowall/2
rY = nY/nColall/2

# Create the X and Y meshes
centerX = collect(range(rX, nX-rX, nRowall));
centerY = collect(range(rY, nY-rY, nColall));

X = collect(range(0.5,nX-0.5,nX)) .* ones(1, nY);
Y = ones(nX) .* collect(range(0.5, nY-0.5,nY))';

sigma = 3
r_x = relative_radius*nX/nRowall/2/sigma
r_y = relative_radius*nY/nColall/2/sigma

out = zeros(ComplexF64, (nX, nY, nRowall, nColall));
# Compute the Gaussian distribution for each BX BY box
for bx in 1:nRowall
    for by in 1:nColall
        x = centerX[bx]
        y = centerY[by]
        mag = exp.(-(((X .- x)/r_x).^2 + ((Y .- y)/r_y).^2) ./ 2)
        pha = angle.((X .-  x) .+ im*(Y .- y)).*0
        # pha = sqrt.((X.-x).^2 .+ (Y.-y).^2)./minimum([nX/nRow/2, nY/nCol/2]).^2 .* π
        # pha = (bx+by) % 2 == 0 ? pha : -pha 
        real = cos.(pha) .* mag
        imag = sin.(pha) .* mag
        out[:, :, bx, by] = (real .+ imag*im)
        # out[:, :, bx, by] = mag #.* exp.(im.*(bx*π+by*π))
    end
end

BX = collect(range(1, nRow*(nBlock-1)+1, nBlock)) .* ones(1, nBlock);
BY = ones(nBlock) .* collect(range(1, nCol*(nBlock-1)+1, nBlock))';
BX = Int.(BX)
BY = Int.(BY)



csm = zeros(ComplexF64, (nX, nY, nCoil));
for i = 1:nRow
    for j = 1:nCol
        for (idx_x, idx_y) in zip(BX.+i.-1, BY.+j.-1)
            csm[:,:, (i-1)*nRow+j] += out[:,:, idx_x, idx_y]
        end
    end
end
norm = sqrt.(sum(abs.(csm) .^ 2, dims=3));
csm = csm./ norm;


fig = plt_images(abs.(csm); dim=3, nRow=nRow, nCol=nCol)
fig = plt_images(angle.(csm); dim=3, nRow=nRow, nCol=nCol, vmin=-π, vmax=π)


# fig = plt_images(abs.(out[:, :, 1, :]); dim=3)
# fig = plt_images(angle.(out); dim=3, nRow=nRowall, nCol=nColall)


# fig, axs = plt.subplots(1, 2)
# ax1 = axs[1]
# ax2 = axs[2]
# for col = 1:nCol
#     ax1.plot(abs.(out[nX÷2, :, (nRow÷2-1)*nCol+col]))
# end
# ax2.scatter(collect(1:1:nCoil), abs.(out)[75, 75, :], s=5)
# ax2.scatter(collect(1:1:nCoil), angle.(out)[75, 75, :], s=5)

