# This file is used to include all the necessary files for the Coil-Sensitivity Map (CSM) generating or loading.
include("Birdcage.jl")
include("binary.jl")
include("Real_32cha.jl")

export csm_Birdcage
export csm_Real_32cha
export csm_Fan_binary
export csm_Rect_binary

include("load_csm.jl")
export load_csm



"""
    factor_a, factor_b = get_factors(num::Int64)
    Get the factors of a number, Where a * b = num, make the difference between a and b as small as possible.

# Arguments
- `num::Int64`: The number whose factors are to be found.

# Returns
- `Tuple{Int64,Int64}`: A tuple containing the two factors of the number.
"""
function get_factors(num::Int64)::Tuple{Int64,Int64}
    num_sqrt = Int64(ceil(sqrt(num)))
    factor_a = num
    factor_b = 1

    for a in range(num_sqrt,1, step=-1)
        if num%a == 0
            factor_a = Int64(a)
            factor_b = Int64(num/a)
            break
        end
    end
    return factor_a, factor_b
end
