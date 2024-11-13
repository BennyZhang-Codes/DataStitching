

"""
    mask = csm_Fan_binary(nX::Int64, nY::Int64, nCoil::Int64; overlap::Real=1, verbose::Bool=false)

# Description
    get fan-shaped csm with arbitrary overlap.

# Arguments
- `nX`: (`::Int64`) 
- `nY`: (`::Int64`) 
- `nCoil`: (`::Int64`) 

# Keywords
- `overlap`: (`::Real`) overlapping factor 
- `verbose`: (`::Real`) print information (default: false)

# Returns
- `mask`: (`::Array{Float64, 3}`) mask (nX, nY, nCoil)

# Examples
```julia-repl
julia> mask = csm_Fan_binary(100, 100, 6; overlap=0.2)
julia> plot_imgs_subplots(mask, 2, 3)
```
"""
function csm_Fan_binary(nX::Int64, nY::Int64, nCoil::Int64; overlap::Real=0.5, verbose::Bool=false)
    if verbose
        @info "fan-shaped binary sensitivity" nX=nX nY=nY nCoil=nCoil overlap=overlap
    end
    m_x = (1:1:nX) .* ones(1, nY)
    m_y = ones(nX) .* (1:1:nY)'

    Δx = vec(m_x) .- (nX/2)
    Δy = vec(m_y) .- (nY/2)

    ϕ = reshape(angle.(Δx .+ Δy*im), (nX, nY))     # split cartesien grids by angle ϕ
    # plot_image(ϕ; zmin=-pi)

    mask = zeros(nX, nY, nCoil)
    for i = 1:nCoil
        angle_rad = collect(-pi:2pi/nCoil:pi)[i]
        m = mask[:,:,i]
        m[abs.(ϕ .- angle_rad) .>= (2pi - pi/nCoil*(1+overlap)) .|| abs.(ϕ .- angle_rad) .<= pi/nCoil*(1+overlap)] .= 1
        mask[:,:,i] = m
    end
    return mask
    # plot_imgs_subplots(mask, 2,2)
end