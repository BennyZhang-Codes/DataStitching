
"""
    mask = csm_Rect_binary(Nx::Int64, Ny::Int64, nCoil::Int64; verbose::Bool=false)

# Description
    get Rect-shaped csm.

# Arguments
- `Nx`: (`::Int64`) 
- `Ny`: (`::Int64`) 
- `nCoil`: (`::Int64`) 

# Keywords
- `verbose`: (`::Real`) print information (default: false)

# Returns
- `mask`: (`::Array{Float64, 3}`) mask (Nx, Ny, nCoil)

# Examples
```julia-repl
julia> mask = csm_Rect_binary(100, 100, 6; verbose=true)
julia> nX = nY = 100
julia> nRow = 3
julia> nCol = 4
julia> mask = csm_Rect_binary(nX, nY, nRow*nCol; nRow=nRow, nCol=nCol, verbose=true)
julia> plt_images(mask; dim=3, nRow=nRow, nCol=nCol)
```
"""
function csm_Rect_binary(
    Nx::Int64, 
    Ny::Int64, 
    nCoil::Int64;     
    nRow=nothing,
    nCol=nothing,
    verbose::Bool=false)

    if nRow === nothing || nCol === nothing
        Npartsx, Npartsy = get_factors(nCoil)
    else
        @assert nRow * nCol == nCoil "nCoil should be equal to nRow * nCol"
        Npartsx = nRow
        Npartsy = nCol
    end

    if verbose
        @info "rectangular binary sensitivity" Nx=Nx Ny=Ny Npartsx=Npartsx Npartsy=Npartsy
    end

    rangex = Int64.(round.(collect(range(0, Nx, Npartsx+1))))
    rangey = Int64.(round.(collect(range(0, Ny, Npartsy+1))))
    mask = zeros(Nx, Ny, Npartsx * Npartsy)
    # [println("$(rangex[i]+1):$(rangex[i+1]), $(rangey[j]+1):$(rangey[j+1]):") for i in eachindex(rangex[1:end-1]), j in eachindex(rangey[1:end-1])]
    for i in eachindex(rangex[1:end-1])
        for j in eachindex(rangey[1:end-1])
            mask[rangex[i]+1:rangex[i+1], rangey[j]+1:rangey[j+1], (i-1)*Npartsy+j] .= 1
        end
    end
    return mask
    # plot_imgs_subplots(mask, Npartsx, Npartsy)
end


function csm_Rect_gaussian(
    Nx::Int64, 
    Ny::Int64, 
    nCoil::Int64;     
    nRow=nothing,
    nCol=nothing,
    verbose::Bool=false)

    if nRow === nothing || nCol === nothing
        Npartsx, Npartsy = get_factors(nCoil)
    else
        @assert nRow * nCol == nCoil "nCoil should be equal to nRow * nCol"
        Npartsx = nRow
        Npartsy = nCol
    end

    if verbose
        @info "rectangular binary sensitivity" Nx=Nx Ny=Ny Npartsx=Npartsx Npartsy=Npartsy
    end

    rangex = Int64.(round.(collect(range(0, Nx, Npartsx+1))))
    rangey = Int64.(round.(collect(range(0, Ny, Npartsy+1))))

    # Define the parameters of the Gaussian distribution
    mean = [0, 0]  # mean of the x and y coordinates
    cov = [[1, 0.5], [0.5, 1]]  # covariance matrix
    m_x = collect(range(-2,2,Nx)) .* ones(1, Ny)
    m_y = ones(Nx) .* collect(range(-2,2,Ny))'
    Z = exp.(-((m_x .- mean[1]).^2 + (m_y .- mean[2]).^2) ./ (2 * cov[1][1]))
    ϕ = exp.(1im.*(Z.*2π.-π))

    mask = zeros(ComplexF64, Nx, Ny, Npartsx * Npartsy)
    # [println("$(rangex[i]+1):$(rangex[i+1]), $(rangey[j]+1):$(rangey[j+1]):") for i in eachindex(rangex[1:end-1]), j in eachindex(rangey[1:end-1])]
    for i in eachindex(rangex[1:end-1])
        for j in eachindex(rangey[1:end-1])
            mask[rangex[i]+1:rangex[i+1], rangey[j]+1:rangey[j+1], (i-1)*Npartsy+j] .= ϕ[rangex[i]+1:rangex[i+1], rangey[j]+1:rangey[j+1]]
        end
    end
    return mask
    # plot_imgs_subplots(mask, Npartsx, Npartsy)
end
