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