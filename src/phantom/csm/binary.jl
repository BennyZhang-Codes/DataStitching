

"""
    mask = csm_Fan_binary(Nx::Int64, Ny::Int64, nCoil::Int64; overlap::Real=1, verbose::Bool=false)

# Description
    get fan-shaped csm with arbitrary overlap.

# Arguments
- `Nx`: (`::Int64`) 
- `Ny`: (`::Int64`) 
- `nCoil`: (`::Int64`) 

# Keywords
- `overlap`: (`::Real`) overlapping factor 
- `verbose`: (`::Real`) print information (default: false)

# Returns
- `mask`: (`::Array{Float64, 3}`) mask (Nx, Ny, nCoil)

# Examples
```julia-repl
julia> mask = csm_Fan_binary(100, 100, 6; overlap=0.2)
julia> plot_imgs_subplots(mask, 2, 3)
```
"""
function csm_Fan_binary(Nx::Int64, Ny::Int64, nCoil::Int64; overlap::Real=0.5, verbose::Bool=false)
    if verbose
        @info "fan-shaped binary sensitivity" Nx=Nx Ny=Ny nCoil=nCoil overlap=overlap
    end
    m_x = (1:1:Nx) .* ones(1, Ny)
    m_y = ones(Nx) .* (1:1:Ny)'

    Δx = vec(m_x) .- (Nx/2)
    Δy = vec(m_y) .- (Ny/2)

    ϕ = reshape(angle.(Δx .+ Δy*im), (Nx, Ny))     # split cartesien grids by angle ϕ
    # plot_image(ϕ; zmin=-pi)

    mask = zeros(Nx, Ny, nCoil)
    for i = 1:nCoil
        angle_rad = collect(-pi:2pi/nCoil:pi)[i]
        m = mask[:,:,i]
        m[abs.(ϕ .- angle_rad) .>= (2pi - pi/nCoil*(1+overlap)) .|| abs.(ϕ .- angle_rad) .<= pi/nCoil*(1+overlap)] .= 1
        mask[:,:,i] = m
    end
    return mask
    # plot_imgs_subplots(mask, 2,2)
end

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


function csm_Rect_gauss_phase(
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


function csm_Ring_binary(Nx::Int64, Ny::Int64, nCoil::Int64; overlap::Real=0.5, verbose::Bool=false)
    if verbose
        @info "Ring-shaped binary sensitivity" Nx=Nx Ny=Ny nCoil=nCoil overlap=overlap
    end
    m_x = (-Nx/2+0.5:1:Nx/2-0.5) .* ones(1, Ny)
    m_y = ones(Nx) .* (-Ny/2+0.5:1:Ny/2-0.5)'
    m_radius = sqrt.(m_x.^2 + (Nx/Ny*m_y).^2)

    Δx = vec(m_x) .- (Nx/2)
    Δy = vec(m_y) .- (Ny/2)

    mask = zeros(Nx, Ny, nCoil)
    rs = collect(0:sqrt(Nx^2+Ny^2)/2/(nCoil-1):sqrt(Nx^2+Ny^2))
    for i = 1:nCoil
        radius_min = rs[i]
        radius_max = rs[i+1]
        m = mask[:,:,i]
        m[m_radius.>= radius_min .&& m_radius.<= radius_max] .= 1
        if i == nCoil
            m[m_radius.>= radius_max] .= 1
        end
        mask[:,:,i] = m
    end
    # plt_images(mask; dim=3)
    return mask
end