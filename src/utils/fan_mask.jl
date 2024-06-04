export get_fan_mask



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
julia> mask = get_fan_mask(100, 100, 6; overlap=1.5)
julia> plot_imgs_subplots(mask, 2, 3)
```
"""
function get_fan_mask(Nx::Int64, Ny::Int64, Nparts::Int64; overlap::Real=1)
    @info "fan mask" Nx=Nx Ny=Ny Nparts=Nparts overlap=overlap
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
        m[abs.(ϕ .- angle_rad) .>= (2pi - pi/Nparts*overlap) .|| abs.(ϕ .- angle_rad) .<= pi/Nparts*overlap] .= 1
        mask[:,:,i] = m
    end
    return mask
    # plot_imgs_subplots(mask, 2,2)
end

