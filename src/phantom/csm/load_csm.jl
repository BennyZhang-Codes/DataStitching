csm_list = [
    :fan, 
    :ring, 
    :rect, 
    :rect_gaussian, 
    :birdcage, 
    :real_32cha, 
    :gaussian_grid, 
    :gaussian_grid_block,
    :gaussian_grid_block_pha,]
export csm_list

function load_csm(
    type::Symbol                    , 
    nX::Int64                       , 
    nY::Int64                       , 
    nCoil::Int64                    ;
    nRow                  = nothing ,
    nCol                  = nothing ,
    nBlock                = 3       ,
    overlap::Real         = 1       ,  # overlap between fan coils, for csm_Fan_binary
    relative_radius::Real = 1.5     ,  # relative radius of the coil, for csm_Birdcage
    use_gpu               = false   ,
    verbose::Bool         = false
)
    @assert type in csm_list "type must be one of the following: $(csm_list)"
    if type == :fan
        csm = csm_Fan_binary(nX, nY, nCoil; overlap=overlap, verbose=verbose)
    elseif type == :ring
        csm = csm_Ring_binary(nX, nY, nCoil; overlap=overlap, verbose=verbose)
    elseif type == :rect
        csm = csm_Rect_binary(nX, nY, nCoil; verbose=verbose, nRow=nRow, nCol=nCol)
    elseif type == :rect_gaussian
        csm = csm_Rect_gaussian(nX, nY, nCoil; verbose=verbose, nRow=nRow, nCol=nCol)
    elseif type == :birdcage
        csm = csm_Birdcage(nX, nY, nCoil; relative_radius=Float64(relative_radius), verbose=verbose)
    elseif type == :real_32cha
        csm = csm_Real_32cha(nX, nY; verbose=verbose)
    elseif type == :gaussian_grid
        csm = csm_Gaussian_grid(nX, nY, nCoil; relative_radius=Float64(relative_radius), nRow=nRow, nCol=nCol, verbose=verbose)
    elseif type == :gaussian_grid_block
        csm = csm_Gaussian_grid_block(nX, nY, nCoil, use_gpu; relative_radius=Float64(relative_radius), nRow=nRow, nCol=nCol, nBlock=nBlock,verbose=verbose)
    elseif type == :gaussian_grid_block_pha
        csm = csm_Gaussian_grid_block_pha(nX, nY, nCoil, use_gpu; relative_radius=Float64(relative_radius), nRow=nRow, nCol=nCol, nBlock=nBlock,verbose=verbose)
    end
    return csm
end


function load_csm(
    type     :: Symbol                      , 
    nX       :: Int64                       , 
    nY       :: Int64                       ,
    nZ       :: Int64                       ,
    nCoil    :: Int64                       ;
    axis     :: String         = "axial"    ,  # orientation
    ss       :: Int64          = 4          ,  # undersample
    location :: Vector{Int64}  = [160, 200] ,  # [slice] for one slice, [start, end] for multiple slices
    verbose  :: Bool           = false      ,
)
    @assert type in csm_list "type must be one of the following: :real_32cha"
    if type == :real_32cha
        csm = csm_Real_32cha(nX, nY, nZ; axis=axis, ss=ss, location=location, verbose=verbose)
    end
    return csm
end