mrddims = ["nCha", "nZ", "nY", "nX", "nAvg", "nSli", "nCon", "nPha", "nRep", "nSet", "nSeg"]
export mrddims

include("get_size.jl")
export get_isize, get_ksize

include("get_data.jl")
export get_kdata