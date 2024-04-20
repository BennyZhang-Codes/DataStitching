export get_fan_mask

function get_fan_mask(Nx::Int64, Ny::Int64, Nparts::Int64)
    m_x = (1:1:Nx) .* ones(1, Ny)
    m_y = ones(Nx) .* (1:1:Ny)'

    Δx = vec(m_x) .- (Nx/2)
    Δy = vec(m_y) .- (Ny/2)

    ϕ = reshape(angle.(Δx .+ Δy*im), (Nx, Ny))     # split cartesien grids by angle ϕ
    # plot_image(ϕ; zmin=-pi)

    mask = zeros(Nx, Ny)
    for i = 1:Nparts
        angle_rad = collect(-pi:2pi/Nparts:pi)[i]
        mask[abs.(ϕ .- angle_rad) .>= (2pi - pi/Nparts) .|| abs.(ϕ .- angle_rad) .<= pi/Nparts] .= i
    end
    return mask
    # plot_image(mask; zmin=minimum(mask), zmax=maximum(mask))
end





