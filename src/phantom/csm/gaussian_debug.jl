nX = 150
nY = 150 
nCoil =25
nRow = 5
nCol = 5
relative_radius = 1
verbose = false

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

# Create the X and Y meshes
X = collect(range(1,nX/rX,nX)) .* ones(1, nY);
Y = ones(nX) .* collect(range(1,nY/rY,nY))';

# Define the mean point coordinates for each BX BY box
meanX = collect(range(1,nX/rX,nRow))
meanY = collect(range(1,nY/rY,nCol))'

# Initialize the Z matrix with zeros
out = zeros((nX, nY, nRow*nCol));
# Define the covariance matrix
cov = [1 0.5; 0.5 1]
# Compute the Gaussian distribution for each BX BY box
for bx in 1:nRow
    for by in 1:nCol
        x = meanX[bx]
        y = meanY[by]
        # Apply the bias term if bx is even
        if by % 2 == 0
            x += nX/rX/nRow/2
        end
        # Compute the Gaussian distribution values for each X and Y coordinate
        out[:, :, (bx-1)*nCol+by] = exp.(-((X .- x).^2 + (Y .- y).^2) ./ (2 * cov[1, 1]))
    end
end

norm = sqrt.(sum(abs.(out) .^ 2, dims=3))
out = out./ norm

fig = plt_images(out; dim=3, nRow=nRow, nCol=nCol)
# fig.savefig("gaussian.png")
# plt_image(sqrt.(sum(out.^2; dims=3))[:,:,1])




