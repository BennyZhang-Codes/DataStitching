"""
    csm_Birdcage(N::Int64, nCoil::Int64, relative_radius::Float64)

# Description
    Computes the sensitivity maps for each coils that are arranged in a birdcage manner.

# Arguments
- `N::Int64`: The number of pixels in the x-direction.
- `nCoil::Int64`: The number of coils in the birdcage.
- `relative_radius::Float64`: The relative radius.

# Returns
   A complex array of size (nX, nY, nCoil) containing the sensitivity maps.
# Examples
```julia-repl
julia> smap = csm_Birdcage(300, 300, 32, 1.5)
julia> plot_imgs_subplots(abs.(smap), 4, 8)
```
"""
function csm_Birdcage(nX::Int64, nY::Int64, nCoil::Int64; relative_radius::Float64=1.5, verbose::Bool=false)
    if verbose
        @info nX=nX nY=nY Npartsx=nCoil
    end

    out = zeros(ComplexF64, nX, nY, nCoil)
    for c=0:nCoil-1
        coilx = relative_radius*cos(c*(2*pi/nCoil))
        coily = relative_radius*sin(c*(2*pi/nCoil))
        coil_phase = -c*(2*pi/nCoil)

        for y=0:nY-1
        y_co = (y - (nY/2))/(nY/2) - coily
        for x=0:nX-1
            x_co = (x - (nX/2))/(nX/2) - coilx
            rr = sqrt(x_co^2+y_co^2)
            phi = atan(x_co, -y_co) + coil_phase
            out[x+1,y+1, c+1] = 1/(rr) * exp(1im*phi)
        end
        end
    end
    norm = sqrt.(sum(abs.(out) .^ 2, dims=3))
    return out./norm
end



# # -*- coding: utf-8 -*-
# """MRI simulation functions.
# """
# import numpy as np

# __all__ = ["birdcage_maps"]


# function birdcage_maps(nX::Int64, nY::Int64, nCoil::Int64, relative_radius::Float64)
#     """Simulates birdcage coil sensitivies.

#     Args:
#         shape (tuple of ints): sensitivity maps shape,
#             can be of length 3, and 4.
#         r (float): relative radius of birdcage.
#         nzz (int): number of coils per ring.
#         dtype (Dtype): data type.

#     Returns:
#         array.
#     """


#     nc, nY, nX
#     c, y, x = np.mgrid[:nc, :nY, :nX]

#     coilx = r * np.cos(c * (2 * np.pi / nc))
#     coily = r * np.sin(c * (2 * np.pi / nc))
#     coil_phs = -c * (2 * np.pi / nc)

#     x_co = (x - nX / 2.0) / (nX / 2.0) - coilx
#     y_co = (y - nY / 2.0) / (nY / 2.0) - coily
#     rr = np.sqrt(x_co**2 + y_co**2)
#     phi = np.arctan2(x_co, -y_co) + coil_phs
#     out = (1.0 / rr) * np.exp(1j * phi)

#     # elif len(shape) == 4:
#     #     nc, nz, nY, nX = shape
#     #     c, z, y, x = np.mgrid[:nc, :nz, :nY, :nX]

#     #     coilx = r * np.cos(c * (2 * np.pi / nzz))
#     #     coily = r * np.sin(c * (2 * np.pi / nzz))
#     #     coilz = np.floor(c / nzz) - 0.5 * (np.ceil(nc / nzz) - 1)
#     #     coil_phs = -(c + np.floor(c / nzz)) * (2 * np.pi / nzz)

#     #     x_co = (x - nX / 2.0) / (nX / 2.0) - coilx
#     #     y_co = (y - nY / 2.0) / (nY / 2.0) - coily
#     #     z_co = (z - nz / 2.0) / (nz / 2.0) - coilz
#     #     rr = (x_co**2 + y_co**2 + z_co**2) ** 0.5
#     #     phi = np.arctan2(x_co, -y_co) + coil_phs
#     #     out = (1 / rr) * np.exp(1j * phi)
#     # else:
#     #     raise ValueError("Can only generate shape with length 3 or 4")

#     rss = sum(abs(out) ** 2, 0) ** 0.5
#     out /= rss

#     return out.astype(dtype)


