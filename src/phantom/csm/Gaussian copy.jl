function csm_Gaussian_grid(
    nX::Int64, 
    nY::Int64, 
    nCoil::Int64;     
    nRow=nothing,
    nCol=nothing,
    relative_radius::Float64=1.,
    verbose::Bool=false)

    if nRow === nothing || nCol === nothing
        nRow, nCol = get_factors(nCoil)
    else
        @assert nRow * nCol == nCoil "nCoil should be equal to nRow * nCol"
    end

    if verbose
        @info "Gaussian grid sensitivity" nX=nX nY=nY nRow=nRow nCol=nCol relative_radius=relative_radius
    end

    rX = nX/nRow/2 * relative_radius
    rY = nY/nCol/2 * relative_radius
    rX_pha = nX/nRow/2 *1
    rY_pha = nY/nCol/2 *1

    # Create the X and Y meshes
    X = collect(range(1,nX/rX,nX)) .* ones(1, nY)
    Y = ones(nX) .* collect(range(1,nY/rY,nY))'
    X_pha = collect(range(1,nX/rX_pha,nX)) .* ones(1, nY);
    Y_pha = ones(nX) .* collect(range(1,nY/rY_pha,nY))';

    # Define the mean point coordinates for each BX BY box
    meanX = collect(range(1,nX/rX,nRow))
    meanY = collect(range(1,nY/rY,nCol))
    meanX_pha = collect(range(1,nX/rX_pha,nRow))
    meanY_pha = collect(range(1,nY/rY_pha,nCol))'

    # Initialize the Z matrix with zeros
    out = zeros(ComplexF64, (nX, nY, nRow*nCol));
    # Define the covariance matrix
    cov = [1 0.5; 0.5 1]
    # Compute the Gaussian distribution for each BX BY box
    for bx in 1:nRow
        for by in 1:nCol
            x = meanX[bx]
            y = meanY[by]
            x_pha = meanX_pha[bx]
            y_pha = meanY_pha[by]
            # Apply the bias term if bx is even
            # if by % 2 == 0
            #     x += nX/rX/nRow/2
            #     x_pha += nX/rX_pha/nRow/2
            # end
            # Compute the Gaussian distribution values for each X and Y coordinate
            mag = exp.(-((X .- x).^2 + (Y .- y).^2) ./ (2 * cov[1, 1]))

            # ϕ = angle.((X .-  x) .+ im*(Y .- y))
            # m = exp.(-((ϕ/2).^2 + (ϕ/2).^2) ./ (2 * cov[1, 1]))
            # co = 1 .- m .* exp.(-((X_pha .- x_pha).^2 + (Y_pha .- y_pha).^2) ./ (2 * cov[1, 1]))
            # pha = ϕ.*co
            # real = cos.(pha) .* mag
            # imag = sin.(pha) .* mag
            # out[:, :, (bx-1)*nCol+by] = (real .+ imag*im)

            out[:, :, (bx-1)*nCol+by] = mag .* exp.(im.*(bx*π+by*π))
        end
    end
    norm = sqrt.(sum(abs.(out) .^ 2, dims=3))
    out = out./ norm
    return out
end




