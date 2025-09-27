b0_list = [
    :real,    
]
export b0_list

"""
    b0map = load_b0map(type::Symbol, nX::Int64, nY::Int64, nZ::Int64; axis="axial", ss=4, location=[160,200], verbose=false)

# Description
    Loads the B0 map of the phantom and resizes it according to the specified dimensions `(nX, nY, nZ)`.  
    You can specify the orientation (`axis`), undersampling factor (`ss`), and slice selection (`location`).  

# Arguments
- `type::Symbol`: type of B0 map to load
- `nX::Int64`: number of voxels of the phantom in the x-direction  
- `nY::Int64`: number of voxels of the phantom in the y-direction  
- `nZ::Int64`: number of voxels of the phantom in the z-direction 
- `axis::String="axial"`: orientation of slices, must be one of `"axial"`, `"coronal"`, `"sagittal"`  
- `ss::Int64=4`: undersampling factor for spatial dimensions  
- `location::Vector{Int64}=[160,200]`:  
    - if length is 1 → index of a single slice  
    - if length is 2 → `[start, end]` slice range  
- `verbose::Bool=false`: whether to print progress messages  

# Returns
- `b0map::Array{Float64, 2}`: The b0map.
"""
function load_b0map(
    type     :: Symbol                      , 
    nX       :: Int64                       , 
    nY       :: Int64                       ,
    nZ       :: Int64                       ;
    axis     :: String         = "axial"    ,  # orientation
    ss       :: Int64          = 4          ,  # undersample
    location :: Vector{Int64}  = [160, 200] ,  # [slice] for one slice, [start, end] for multiple slices
    verbose  :: Bool           = false      ,
)
    @assert type in b0_list "type must be one of the following: :real"
    if type == :real
        b0 = b0_real(nX, nY, nZ; axis=axis, ss=ss, location=location, verbose=verbose)
    end
    return b0
end
