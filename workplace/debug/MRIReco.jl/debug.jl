F = FFTOp(ComplexF64, shape=size(img), shift=true)
# apply operator
k = F * vec(img)            
# apply the adjoint operator
image = adjoint(F) * k     

K = reshape(k, size(img)...)  # reshape to 2D
plot_image(abs.(K).^0.02)
image = reshape(image, size(img)...)  # reshape to 2D
plot_image(abs.(image))
