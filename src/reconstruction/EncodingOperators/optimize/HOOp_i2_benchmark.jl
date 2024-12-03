using BenchmarkTools
solver = "cgnr"; regularization = "L2"; λ = 1.e-4; iter=10;
recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
recParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, csm_nCoil));

numContr, numChan = MRIReco.numContrasts(acqData), MRIReco.numChannels(acqData);
reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, recParams);
senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan);
smaps = senseMaps[:,:,1,:];
S = SensitivityOp(reshape(ComplexF64.(smaps),:,numChan),1)

HOOp = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); 
                        Nblocks=10, fieldmap=Matrix(B0map), grid=1, use_gpu=true, verbose=true);
Op   = DiagOp(HOOp, nCha) ∘ S 
Opi2 = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); 
                        Nblocks=10, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=Matrix(B0map), grid=1, use_gpu=true, verbose=true);
                        
W = WeightingOp(Complex{T}; weights=weights[1], rep=numChan)
E1 = ∘(W, Op, isWeighting=true)
E2 = ∘(W, Opi2, isWeighting=false)
EH1 = normalOperator(E1)
EH2 = normalOperator(E2)

recoParams = merge(defaultRecoParams(), recParams)
solver1 = createLinearSolver("cgnr", E1; AᴴA=EH1, recoParams...)
solver2 = createLinearSolver("cgnr", E2; AᴴA=EH2, recoParams...)

kdata = multiCoilData(acqData, 1, 1, rep=1) .* repeat(weights[1], numChan)

@benchmark I1 = solve(solver1, kdata)
plt_image(rotl90(abs.(reshape(I1, (Nx, Ny)))))


@benchmark solve(solver2, kdata)
plt_image(rotl90(abs.(reshape(I2, (Nx, Ny)))))

Op = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); 
                        Nblocks=5, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=Matrix(B0map), grid=1, use_gpu=true, verbose=true);


@time x = recon_HOOp_i2(acqData, Op, recParams)
plt_image(rotl90(abs.(x)))

x = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (1., 1.); location=location, ss=simtype.ss);
plt_image(rotl90(x))








# debug
reconsize    = (Nx, Ny)
tr_nominal   = tr_nominal
tr_measured  = tr_dfc_stitched
sim_method   = BlochHighOrder("000")
fieldmap     = Matrix(B0map)
csm          = Complex{T}.(sensitivity[:,:,1,:])
Nblocks      = 10
use_gpu      = true
verbose      = true
Δx           = 1e-3
Δy           = 1e-3
grid         = 1

# init of HighOrderOp
nodes_measured = T.(kspaceNodes(tr_measured))
nodes_nominal = T.(kspaceNodes(tr_nominal))
times = T.(readoutTimes(tr_measured))
@assert size(nodes_measured,1) == 9 "nodes for measured must have 9 rows"
@assert size(nodes_nominal,1) == 3 "nodes for nominal must have 3 rows"
@assert size(fieldmap) == reconsize "fieldmap must have same size as reconsize"
@assert size(csm)[1:2] == reconsize "fieldmap must have same size as reconsize"

nX, nY = reconsize
nCha = size(csm, 3)
nSam = size(nodes_measured,2)
nRow = nSam * nCha
nCol = nVox = prod(reconsize)

csm = reshape(csm, nX*nY, nCha)  # [nX*nY, nCha]
fieldmap = vec(fieldmap)

k = nSam # number of nodes
Nblocks = Nblocks > k ? k : Nblocks # Nblocks must be <= k

n = k÷Nblocks # number of nodes per block
parts = [n for i=1:Nblocks] # number of nodes per block
parts = [1+n*(i-1):n*i for i=1:Nblocks]
if k%Nblocks!= 0
    push!(parts, n*Nblocks+1:k)
end


x, y = 1:nX, 1:nY
if grid == 1      # x up->down, y left->right
    x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- nX/2 .- 1, y .- nY/2 .- 1
elseif grid == 2
    x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- (nX+1)/2, y .- (nY+1)/2
elseif grid == 3  # x left->right, y down->up
    x, y, z = vec(x*0.0 .+ y'), vec(x[end:-1:1] .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- (nX+1)/2, y .- (nY+1)/2
elseif grid == 4  # x left->right, y up->down
    x, y, z = vec(x*0.0 .+ y'), vec(x .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- (nX+1)/2, y .- (nY+1)/2
elseif grid == 5  # x right->left, y up->down
    x, y, z = vec(x*0.0 .+ y[end:-1:1]'), vec(x .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- (nX+1)/2, y .- (nY+1)/2
end
# print(x)
x, y = x * Δx, y * Δy 

# input data
x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (1., 1.); location=location, ss=5);
xm = vec(x_ref)

Nblocks = 10
function prod_HighOrderOp_i3(
    xm::AbstractVector{T}, 
    x::Vector{Float64}, 
    y::Vector{Float64}, 
    z::Vector{Float64}, 
    nVox::Int64, 
    nSam::Int64, 
    nCha::Int64,
    nodes_measured::Matrix{Float64}, 
    nodes_nominal::Matrix{Float64},
    times::Vector{Float64}, 
    fieldmap::Vector{Float64},
    csm::Array{ComplexF64, 2};
    sim_method::BlochHighOrder=BlochHighOrder("111"), 
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nSam], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where T<:Union{Real,Complex}
    xm = Vector(xm)
    if verbose
        @info "HighOrderOp_i2 prod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    out = zeros(ComplexF64, nSam, nCha)
    x0 = ones(Float64, nVox)
    if use_gpu
        out = out |> gpu
        nodes_measured = nodes_measured |> gpu
        nodes_nominal = nodes_nominal |> gpu
        xm = xm |> gpu
        x0 = x0 |> gpu
        x = x |> gpu
        y = y |> gpu
        z = z |> gpu
        times = times |> gpu
        fieldmap = fieldmap |> gpu
        csm = csm |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        @inbounds begin
            h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(nodes_measured[1,p]), @view(nodes_measured[2,p]), @view(nodes_measured[3,p]), @view(nodes_measured[4,p]),
                                        @view(nodes_measured[5,p]), @view(nodes_measured[6,p]), @view(nodes_measured[7,p]), @view(nodes_measured[8,p]), @view(nodes_measured[9,p])
            hx, hy, hz = @view(nodes_nominal[1,p]), @view(nodes_nominal[2,p]), @view(nodes_nominal[3,p])
            ϕ0 = sim_method.ho0 ? h0 .* x0' : 0
            ϕ1 = sim_method.ho1 ? (h1 .* x') .+ (h2 .* y') .+ (h3 .* z') : (hx .* x') .+ (hy .* y') .+ (hz .* z')
            ϕ2 = sim_method.ho2 ? h4 .* (x .* y)' .+ h5 .* (z .* y)' .+ h6 .* (3z.^2-(x.^2 .+ y.^2 .+ z.^2))' .+
                    h7 .* (x .* z)' .+ h8 .* (x.^2 .- y.^2)' : 0
            ϕB0 = times[p] .* fieldmap'
            ϕ = ϕ0 .+ ϕ1 .+ ϕ2 .+ ϕB0
            e = exp.(-2*1im*pi*ϕ)
            out[p, :] =  e * (xm .* csm)
        end
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end

@benchmark prod_HighOrderOp_i2(xm, x, y, z, nVox, nSam, nCha, nodes_measured, nodes_nominal, times, fieldmap, csm;
sim_method, Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose)

using KomaHighOrder: prod_HighOrderOp
@benchmark prod_HighOrderOp(xm, x, y, z, nodes_measured, nodes_nominal, times, fieldmap;
                                        sim_method, Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose)



@btime begin
    if verbose
        @info "HighOrderOp_i2 prod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    out = zeros(ComplexF64, nSam, nCha)
    x0 = ones(Float64, nVox)
    if use_gpu
        out = out |> gpu
        nodes_measured = nodes_measured |> gpu
        nodes_nominal = nodes_nominal |> gpu
        xm = xm |> gpu
        x0 = x0 |> gpu
        x = x |> gpu
        y = y |> gpu
        z = z |> gpu
        times = times |> gpu
        fieldmap = fieldmap |> gpu
        csm = csm |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(nodes_measured[1,p]), @view(nodes_measured[2,p]), @view(nodes_measured[3,p]), @view(nodes_measured[4,p]),
                                    @view(nodes_measured[5,p]), @view(nodes_measured[6,p]), @view(nodes_measured[7,p]), @view(nodes_measured[8,p]), @view(nodes_measured[9,p])
        hx, hy, hz = @view(nodes_nominal[1,p]), @view(nodes_nominal[2,p]), @view(nodes_nominal[3,p])
        ϕ0 = sim_method.ho0 ? h0 .* x0' : 0
        ϕ1 = sim_method.ho1 ? (h1 .* x') .+ (h2 .* y') .+ (h3 .* z') : (hx .* x') .+ (hy .* y') .+ (hz .* z')
        ϕ2 = sim_method.ho2 ? h4 .* (x .* y)' .+ h5 .* (z .* y)' .+ h6 .* (3z.^2-(x.^2 .+ y.^2 .+ z.^2))' .+
                h7 .* (x .* z)' .+ h8 .* (x.^2 .- y.^2)' : 0
        ϕB0 = times[p] .* fieldmap'
        ϕ = ϕ0 .+ ϕ1 .+ ϕ2 .+ ϕB0
        e = exp.(-2*1im*pi*ϕ)

        # @time out[p, 1] =  e * xm 
        # @time out[p, :] =  e * (xm .* csm)
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end
end



xm = vec(x_ref)
@time begin
    if verbose
        @info "HighOrderOp_i2 prod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    # xm = repeat(xm, 1, nCha)
    out = zeros(ComplexF64, nSam, nCha)
    x0 = ones(Float64, nVox)
    if use_gpu
        out = out |> gpu
        nodes_measured = nodes_measured |> gpu
        nodes_nominal = nodes_nominal |> gpu
        xm = xm |> gpu
        x0 = x0 |> gpu
        x = x |> gpu
        y = y |> gpu
        z = z |> gpu
        times = times |> gpu
        fieldmap = fieldmap |> gpu
        csm = csm |> gpu
    end
    # progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(nodes_measured[1,p]), @view(nodes_measured[2,p]), @view(nodes_measured[3,p]), @view(nodes_measured[4,p]),
                                       @view(nodes_measured[5,p]), @view(nodes_measured[6,p]), @view(nodes_measured[7,p]), @view(nodes_measured[8,p]), @view(nodes_measured[9,p])
        hx, hy, hz = @view(nodes_nominal[1,p]), @view(nodes_nominal[2,p]), @view(nodes_nominal[3,p])
        ϕ0 = sim_method.ho0 ? h0 .* x0' : 0
        ϕ1 = sim_method.ho1 ? (h1 .* x') .+ (h2 .* y') .+ (h3 .* z') : (hx .* x') .+ (hy .* y') .+ (hz .* z')
        ϕ2 = sim_method.ho2 ? h4 .* (x .* y)' .+ h5 .* (z .* y)' .+ h6 .* (3z.^2-(x.^2 .+ y.^2 .+ z.^2))' .+
                h7 .* (x .* z)' .+ h8 .* (x.^2 .- y.^2)' : 0
        ϕB0 = times[p] .* fieldmap'
        ϕ = ϕ0 .+ ϕ1 .+ ϕ2 .+ ϕB0
        e = exp.(-2*1im*pi*ϕ)
        # out[p, :] =  e * (xm .* csm)
        # if verbose
        #     next!(progress_bar, showvalues=[(:Nblocks, block)])
        # end
    end
end


