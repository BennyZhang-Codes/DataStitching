use_gpu = true
Nblocks = 20
shape = (300, 300)

nodes = Float64.(kspaceNodes(tr))
times = readoutTimes(tr)
nrow = size(nodes,2)
ncol = prod(shape)
k = size(nodes,2) # number of nodes
Nblocks = Nblocks > k ? k : Nblocks # Nblocks must be <= k

n = k÷Nblocks # number of nodes per block
parts = [n for i=1:Nblocks] # number of nodes per block
parts = [1+n*(i-1):n*i for i=1:Nblocks]
if k%Nblocks!= 0
    push!(parts, n*Nblocks+1:k)
end

Nx, Ny = shape
x, y = 1:Nx, 1:Ny
x, y = vec(x .+ y'*0), vec(x*0 .+ y') #grid points
x, y = x .- Nx/2 .- 1, y .- Ny/2 .- 1



#  prod
out = zeros(ComplexF64,size(nodes,2))
xm = ones(Float64, prod(shape))
if use_gpu
    out = out |> gpu
    nodes = nodes |> gpu
    xm = xm |> gpu
    x = x |> gpu
    y = y |> gpu
end

@time for (block, p) = enumerate(parts)
    kx, ky = @view(nodes[1,p]), @view(nodes[2,p])
    ϕ = (kx .* x') .+ (ky .* y')
    e = exp.(-2*1im*pi*ϕ)
    out[p] =  e * xm
end
@time e1 = exp.(-2*1im*pi*ϕ);
@time e2 = cos.(-2pi.*ϕ) .+ 1im*sin.(-2pi.*ϕ);

@time ϕ = (kx .* x') .+ (ky .* y');
@time e = exp.(-2*1im*pi*ϕ);
@time out[p] =  e * xm;

@time for (block, p) = enumerate(parts)
    kx, ky = @view(nodes[1,p]), @view(nodes[2,p])
end

@time for (block, p) = enumerate(parts)
    kx, ky = @view(nodes[1,p]), @view(nodes[2,p])
end



##############################
#  expϕ
##############################
zeros(ComplexF64, (size(nodes,2), size(x, 1)))




