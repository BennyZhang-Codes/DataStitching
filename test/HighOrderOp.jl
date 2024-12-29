using Test
using KomaHighOrder
using PyPlot
using MRISimulation


# setting some parameters
T = Float32
nX = nY = 64; nZ = 1
Δx = Δy = Δz = T.(1. * 1e-3)   # [m]
nCha = 8;
nSample = nX * nY;
kmax_x = kmax_y = kmax_z = 1/Δx;  # [m^-1]
kx = collect(range(-0.5, stop=0.5, length=nX+1))[1:nX] * kmax_x;
ky = collect(range(-0.5, stop=0.5, length=nY+1))[1:nY] * kmax_y;
kx = kx .* ones(nY)';
ky = ky' .* ones(nX);

kspha = zeros(9, nSample);
kspha[2, :] = vec(kx);
kspha[3, :] = vec(ky);

grid = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=false, reverse_y=true);
csm = csm_Birdcage(nX, nY, nCha);

b0 = quadraticFieldmap(nX, nY, 50.)[:,:,1];

t = collect(range(-1e-6, stop=1e-6, length=nX+1))[1:nX] * nX
datatime = vec(t .* ones(nY)')

plt_images(  abs.(csm); dim=3, nRow=2, nCol=4)
plt_images(angle.(csm); dim=3, nRow=2, nCol=4)
plt_B0map(b0)

weight = SampleDensity(kspha[2:3, :], (nX, nY));

# init HighOrderOpv2_i2
Nblocks = 10
use_gpu = true 
verbose = true
BHO     = BlochHighOrder("111")
HOOp = HighOrderOp(grid, T.(kspha), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);


# prod & ctprod
x = ones(nX*nY);          # x are all ones
y = HOOp * x;
x1 = HOOp' * y;
x = reshape(x, nX, nY);
x1 = reshape(x1, nX, nY);

# judge if x ≈ x1, with digits setting the accuracy
@test x ≈ round.(abs.(x1), digits=2)

fig = plt.Figure()
subplot(1,3,1);plt.imshow(abs.(x1), cmap="gray");plt.title("x1");plt.axis("off")
subplot(1,3,2);plt.imshow(abs.(x-x1), cmap="gray");plt.title("abs(x-x1)");plt.axis("off")
subplot(1,3,3);plt.imshow(angle.(x1), cmap="gray");plt.title("angle(x1)");plt.axis("off")


W = WeightingOp(Complex{T}; weights=Complex{T}.(weight), rep=nCha)
E = ∘(W, HOOp, isWeighting=false)

x = ones(nX*nY);          # x are all ones
y_HOOp = HOOp  * x;
y_E    = E     * x;
x_E    = E'    * y_E;
x_HOOp = HOOp' * y_HOOp;    

weight1 = y_E ./ y_HOOp;
@test round.(abs.(weight), digits=4) ≈ round.(abs.(weight1[1:nSample]), digits=4)

