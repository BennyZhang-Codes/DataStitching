"""
    csm_Birdcage(N::Int64, ncoils::Int64, relative_radius::Float64)

Computes the sensitivity maps for each coils that are arranged
in a birdcage manner.
# out
   A complex array of size (Nx, Ny, ncoils) containing the sensitivity maps.
# Examples
```julia-repl
julia> smap = csm_Birdcage(300, 300, 32, 1.5)
julia> plot_imgs_subplots(abs.(smap), 4, 8)
```
"""
function csm_Birdcage(Nx::Int64, Ny::Int64, ncoils::Int64; relative_radius::Float64=1.5, verbose::Bool=false)
    if verbose
        @info Nx=Nx Ny=Ny Npartsx=ncoils
    end

    out = zeros(ComplexF64, Nx, Ny, ncoils)
    for c=0:ncoils-1
        coilx = relative_radius*cos(c*(2*pi/ncoils))
        coily = relative_radius*sin(c*(2*pi/ncoils))
        coil_phase = -c*(2*pi/ncoils)

        for y=0:Ny-1
        y_co = (y - (Ny/2))/(Ny/2) - coily
        for x=0:Nx-1
            x_co = (x - (Nx/2))/(Nx/2) - coilx
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


# function birdcage_maps(Nx::Int64, Ny::Int64, ncoils::Int64, relative_radius::Float64)
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


#     nc, ny, nx
#     c, y, x = np.mgrid[:nc, :ny, :nx]

#     coilx = r * np.cos(c * (2 * np.pi / nc))
#     coily = r * np.sin(c * (2 * np.pi / nc))
#     coil_phs = -c * (2 * np.pi / nc)

#     x_co = (x - nx / 2.0) / (nx / 2.0) - coilx
#     y_co = (y - ny / 2.0) / (ny / 2.0) - coily
#     rr = np.sqrt(x_co**2 + y_co**2)
#     phi = np.arctan2(x_co, -y_co) + coil_phs
#     out = (1.0 / rr) * np.exp(1j * phi)

#     # elif len(shape) == 4:
#     #     nc, nz, ny, nx = shape
#     #     c, z, y, x = np.mgrid[:nc, :nz, :ny, :nx]

#     #     coilx = r * np.cos(c * (2 * np.pi / nzz))
#     #     coily = r * np.sin(c * (2 * np.pi / nzz))
#     #     coilz = np.floor(c / nzz) - 0.5 * (np.ceil(nc / nzz) - 1)
#     #     coil_phs = -(c + np.floor(c / nzz)) * (2 * np.pi / nzz)

#     #     x_co = (x - nx / 2.0) / (nx / 2.0) - coilx
#     #     y_co = (y - ny / 2.0) / (ny / 2.0) - coily
#     #     z_co = (z - nz / 2.0) / (nz / 2.0) - coilz
#     #     rr = (x_co**2 + y_co**2 + z_co**2) ** 0.5
#     #     phi = np.arctan2(x_co, -y_co) + coil_phs
#     #     out = (1 / rr) * np.exp(1j * phi)
#     # else:
#     #     raise ValueError("Can only generate shape with length 3 or 4")

#     rss = sum(abs(out) ** 2, 0) ** 0.5
#     out /= rss

#     return out.astype(dtype)


