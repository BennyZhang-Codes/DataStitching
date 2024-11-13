function csm_Ring_binary(nX::Int64, nY::Int64, nCoil::Int64; overlap::Real=0.5, verbose::Bool=false)
    if verbose
        @info "Ring-shaped binary sensitivity" nX=nX nY=nY nCoil=nCoil overlap=overlap
    end
    m_x = (-nX/2+0.5:1:nX/2-0.5) .* ones(1, nY)
    m_y = ones(nX) .* (-nY/2+0.5:1:nY/2-0.5)'
    m_radius = sqrt.(m_x.^2 + (nX/nY*m_y).^2)

    Δx = vec(m_x) .- (nX/2)
    Δy = vec(m_y) .- (nY/2)

    mask = zeros(nX, nY, nCoil)
    rs = collect(0:sqrt(nX^2+nY^2)/2/(nCoil-1):sqrt(nX^2+nY^2))
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