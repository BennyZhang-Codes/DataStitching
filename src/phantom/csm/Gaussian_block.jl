function csm_Gaussian_grid_block(
    nX::Int64, 
    nY::Int64, 
    nCoil::Int64;     
    nRow=nothing,
    nCol=nothing,
    nBlock::Int64=3, 
    relative_radius::Real=1.,
    verbose::Bool=false)

    nRowall  = nRow * nBlock
    nColall  = nCol * nBlock

    if nRow === nothing || nCol === nothing
        nRow, nCol = get_factors(nCoil)
    else
        @assert nRow * nCol == nCoil "nCoil should be equal to nRow * nCol"
    end

    if verbose
        @info "Gaussian grid block sensitivity" nX=nX nY=nY nRow=nRow nCol=nCol nBlock=nBlock nCoil=nCoil nCoil_all=nRowall*nColall relative_radius=relative_radius 
    end

    rX = nX/nRowall/2
    rY = nY/nColall/2
    

    sigma = 3
    r_x = relative_radius*nX/nRowall/2/sigma
    r_y = relative_radius*nY/nColall/2/sigma
    # Create the X and Y meshes
    centerX = collect(range(rX-nX/nBlock, nX*(nBlock+1)/nBlock-rX, nRow*(nBlock+2)));
    centerY = collect(range(rY-nY/nBlock, nY*(nBlock+1)/nBlock-rY, nCol*(nBlock+2)));
    
    X = collect(range(0.5,nX-0.5,nX)) .* ones(1, nY);
    Y = ones(nX) .* collect(range(0.5, nY-0.5,nY))';
    BX = collect(range(1, nRow*(nBlock+1)+1, nBlock+2)) .* ones(1, nBlock+2);
    BY = ones(nBlock+2) .* collect(range(1, nCol*(nBlock+1)+1, nBlock+2))';
    
    csm = zeros(ComplexF64, (nX, nY, nCoil));
    # Compute the Gaussian distribution for each BX BY box
    for row in 1:nRow
        for col in 1:nCol
            for i in 1:nBlock+2
                for j in 1:nBlock+2
                    x = centerX[Int64(BX[i,j])+row-1]
                    y = centerY[Int64(BY[i,j])+col-1]
                    mag = exp.(-(((X .- x)/r_x).^2 + ((Y .- y)/r_y).^2) ./ 2)
                    # pha = angle.((X .-  x) .+ im*(Y .- y)).*0
                    # real = cos.(pha) .* mag
                    # imag = sin.(pha) .* mag
                    # out[:, :, bx, by] = (real .+ imag*im)
                    csm[:, :, (row-1)*nCol+col] += mag.+0im
                end
            end
        end
    end
    norm = sqrt.(sum(abs.(csm) .^ 2, dims=3));
    csm = csm./ norm;
    return csm
end

function csm_Gaussian_grid_block0(
    nX::Int64, 
    nY::Int64, 
    nCoil::Int64;     
    nRow=nothing,
    nCol=nothing,
    nBlock::Int64=3, 
    relative_radius::Real=1.,
    verbose::Bool=false)

    nRowall  = nRow * nBlock
    nColall  = nCol * nBlock

    if nRow === nothing || nCol === nothing
        nRow, nCol = get_factors(nCoil)
    else
        @assert nRow * nCol == nCoil "nCoil should be equal to nRow * nCol"
    end

    if verbose
        @info "Gaussian grid block sensitivity" nX=nX nY=nY nRow=nRow nCol=nCol nBlock=nBlock nCoil=nCoil nCoil_all=nRowall*nColall relative_radius=relative_radius 
    end

    rX = nX/nRowall/2
    rY = nY/nColall/2
    
    # Create the X and Y meshes
    centerX = collect(range(rX, nX-rX, nRowall));
    centerY = collect(range(rY, nY-rY, nColall));
    
    X = collect(range(0.5,nX-0.5,nX)) .* ones(1, nY);
    Y = ones(nX) .* collect(range(0.5, nY-0.5,nY))';
    
    sigma = 3
    r_x = relative_radius*nX/nRowall/2/sigma
    r_y = relative_radius*nY/nColall/2/sigma
    
    out = zeros(ComplexF64, (nX, nY, nRowall, nColall));
    # Compute the Gaussian distribution for each BX BY box
    for bx in 1:nRowall
        for by in 1:nColall
            x = centerX[bx]
            y = centerY[by]
            mag = exp.(-(((X .- x)/r_x).^2 + ((Y .- y)/r_y).^2) ./ 2)
            pha = angle.((X .-  x) .+ im*(Y .- y)).*0
            # pha = sqrt.((X.-x).^2 .+ (Y.-y).^2)./minimum([nX/nRow/2, nY/nCol/2]).^2 .* π
            # pha = (bx+by) % 2 == 0 ? pha : -pha 
            real = cos.(pha) .* mag
            imag = sin.(pha) .* mag
            out[:, :, bx, by] = (real .+ imag*im)
            # out[:, :, bx, by] = mag #.* exp.(im.*(bx*π+by*π))
        end
    end
    
    BX = collect(range(1, nRow*(nBlock-1)+1, nBlock)) .* ones(1, nBlock);
    BY = ones(nBlock) .* collect(range(1, nCol*(nBlock-1)+1, nBlock))';
    BX = Int.(BX);
    BY = Int.(BY);

    csm = zeros(ComplexF64, (nX, nY, nCoil));
    for i = 1:nRow
        for j = 1:nCol
            for (idx_x, idx_y) in zip(BX.+i.-1, BY.+j.-1)
                csm[:,:, (i-1)*nRow+j] += out[:,:, idx_x, idx_y]
            end
        end
    end
    norm = sqrt.(sum(abs.(csm) .^ 2, dims=3));
    csm = csm./ norm;
    return csm
end
