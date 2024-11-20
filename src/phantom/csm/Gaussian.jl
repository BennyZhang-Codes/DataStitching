function csm_Gaussian_grid(
    nX::Int64, 
    nY::Int64, 
    nCoil::Int64;     
    nRow=nothing,
    nCol=nothing,
    relative_radius::Real=1.,
    verbose::Bool=false)

    if nRow === nothing || nCol === nothing
        nRow, nCol = get_factors(nCoil)
    else
        @assert nRow * nCol == nCoil "nCoil should be equal to nRow * nCol"
    end

    if verbose
        @info "Gaussian grid sensitivity" nX=nX nY=nY nRow=nRow nCol=nCol relative_radius=relative_radius
    end


    rX = nX/nRow/2
    rY = nY/nCol/2
    
    # Create the X and Y meshes
    centerX = collect(range(rX, nX-rX, nRow));
    centerY = collect(range(rY, nY-rY, nCol));
    
    X = collect(range(0.5,nX-0.5,nX)) .* ones(1, nY);
    Y = ones(nX) .* collect(range(0.5, nY-0.5,nY))';
    
    sigma = 3
    r_x = relative_radius*nX/nRow/2/sigma
    r_y = relative_radius*nY/nCol/2/sigma
    
    out = zeros(ComplexF64, (nX, nY, nCoil));
    # Compute the Gaussian distribution for each BX BY box
    for bx in 1:nRow
        for by in 1:nCol
            x = centerX[bx]
            y = centerY[by]
    
            mag = exp.(-(((X .- x)/r_x).^2 + ((Y .- y)/r_y).^2) ./ 2)
            pha = angle.((X .-  x) .+ im*(Y .- y))
            #pha = sqrt.((X.-x).^2 .+ (Y.-y).^2)./minimum([nX/nRow/2, nY/nCol/2]).^2 .* π
    
            real = cos.(pha) .* mag
            imag = sin.(pha) .* mag
            out[:, :, (bx-1)*nCol+by] = (real .+ imag*im)
            # out[:, :, (bx-1)*nCol+by] = mag .* exp.(im.*(bx*π+by*π))
        end
    end
    norm = sqrt.(sum(abs.(out) .^ 2, dims=3));
    out = out./ norm;    
    return out
end
