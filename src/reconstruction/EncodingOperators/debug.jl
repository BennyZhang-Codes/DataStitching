Δx = Δy = 1e-3

Nx = Ny = 3
x, y = 1:Nx, 1:Ny
x, y = x .+ y'*0.0, x*0.0 .+ y'
# x, y = x .- Nx/2 .- 1, y .- Ny/2 .- 1
x, y = x .- (Nx+1)/2, y .- (Ny+1)/2
x, y = x * Δx, y * Δy 




a = [1 2 3;4 5 6;7 8 9]
a[get_center_range(3, 2), get_center_range(3, 2)]