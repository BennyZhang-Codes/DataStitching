using KomaHighOrder
using PyPlot, LinearAlgebra

nX    = 500
nY    = 500
nRow = 5
nCol = 5
nCoil = nRow*nCol
nBlock = 3

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

sigma = 3
r_x = relative_radius*nX/nRowall/2/sigma
r_y = relative_radius*nY/nColall/2/sigma
# Create the X and Y meshes
centerX = collect(range(rX-nX/nBlock, nX*(nBlock+1)/nBlock-rX, nRow*(nBlock+2)));
centerY = collect(range(rY-nY/nBlock, nY*(nBlock+1)/nBlock-rY, nCol*(nBlock+2)));

X = collect(range(0.5,nX-0.5,nX)) .* ones(1, nY);
Y = ones(nX) .* collect(range(0.5, nY-0.5,nY))';
BX = collect(range(1, nRow*(nBlock+1)+1, nBlock+2)) .* ones(1, nBlock+2);
BY = ones(nBlock+2) .* collect(range(1, nCol*(nBlock+1)+1, nBlock+2))';

csm = zeros(ComplexF64, (nX, nY, nCoil));
# Compute the Gaussian distribution for each BX BY box
@time for row in 1:nRow, col in 1:nCol, i in 1:nBlock+2, j in 1:nBlock+2
    x = centerX[Int64(BX[i,j])+row-1]
    y = centerY[Int64(BY[i,j])+col-1]
    mag = exp.(-(((X .- x)/r_x).^2 + ((Y .- y)/r_y).^2) ./ 2)
    # pha = angle.((X .-  x) .+ im*(Y .- y)).*0
    # real = cos.(pha) .* mag
    # imag = sin.(pha) .* mag
    # out[:, :, bx, by] = (real .+ imag*im)
    csm[:, :, (row-1)*nCol+col] += mag.+0im
end

# Nthreads=Threads.nthreads()

# parts = kfoldperm(length(obj), Nthreads)
# dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times

# ThreadsX.foreach(enumerate(parts)) do (i, p)
#     run_spin_precession!(@view(obj[p]), hoseqd, @view(sig[dims...,i]), @view(Xt[p]), sim_method)
# end


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

