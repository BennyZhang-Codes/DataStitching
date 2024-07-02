

"""
    mask = get_fan_mask(Nx::Int64, Ny::Int64, Nparts::Int64; overlap::Real=1)

# Description
    get fan-shaped mask with arbitrary overlap.

# Arguments
- `Nx`: (`::Int64`) 
- `Ny`: (`::Int64`) 
- `Nparts`: (`::Int64`) 

# Keywords
- `overlap`: (`::Real`) overlapping factor 

# Returns
- `mask`: (`::Array{Float64, 3}`) mask (Nx, Ny, Nparts)

# Examples
```julia-repl
julia> mask = get_fan_mask(100, 100, 6; overlap=0.2)
julia> plot_imgs_subplots(mask, 2, 3)
```
"""
function get_fan_mask(Nx::Int64, Ny::Int64, Nparts::Int64; overlap::Real=0.5, verbose::Bool=false)
    if verbose
        @info "fan-shaped binary sensitivity" Nx=Nx Ny=Ny Nparts=Nparts overlap=overlap
    end
    m_x = (1:1:Nx) .* ones(1, Ny)
    m_y = ones(Nx) .* (1:1:Ny)'

    Δx = vec(m_x) .- (Nx/2)
    Δy = vec(m_y) .- (Ny/2)

    ϕ = reshape(angle.(Δx .+ Δy*im), (Nx, Ny))     # split cartesien grids by angle ϕ
    # plot_image(ϕ; zmin=-pi)

    mask = zeros(Nx, Ny, Nparts)
    for i = 1:Nparts
        angle_rad = collect(-pi:2pi/Nparts:pi)[i]
        m = mask[:,:,i]
        m[abs.(ϕ .- angle_rad) .>= (2pi - pi/Nparts*(1+overlap)) .|| abs.(ϕ .- angle_rad) .<= pi/Nparts*(1+overlap)] .= 1
        mask[:,:,i] = m
    end
    return mask
    # plot_imgs_subplots(mask, 2,2)
end

function get_rect_mask(Nx::Int64, Ny::Int64, Npartsx::Int64, Npartsy::Int64, verbose::Bool=false)
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

# Nx = 150
# Ny = 150
# Npartsx = 5
# Npartsy = 6
# rangex = Int64.(round.(collect(range(0, Nx, Npartsx+1))))
# rangey = Int64.(round.(collect(range(0, Ny, Npartsy+1))))
# mask = zeros(Nx, Ny, Npartsx * Npartsy)
# [println("$(rangex[i]+1):$(rangex[i+1]), $(rangey[j]+1):$(rangey[j+1]):") for i in eachindex(rangex[1:end-1]), j in eachindex(rangey[1:end-1])]
# for i in eachindex(rangex[1:end-1])
#     for j in eachindex(rangey[1:end-1])
#         mask[rangex[i]+1:rangex[i+1], rangey[j]+1:rangey[j+1], (i-1)*Npartsy+j] .= 1
#     end
# end
# plot_imgs_subplots(mask, Npartsx, Npartsy)
