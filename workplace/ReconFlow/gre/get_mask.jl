function get_mask(image; threshold=0)
    image = abs.(image);
    mask = zeros(size(image));
    mask[image.>threshold] .= 1;
    mask = isone.(mask);
    return mask
end
