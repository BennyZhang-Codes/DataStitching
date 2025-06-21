using MRIFieldmaps: b0scale, b0init, b0map
export get_B0map, get_mask

function get_B0map(ydata, TE, smap, mask; l2b=0.002, niter=40, precon=:diag)
    # MRIFieldmaps.jl - 3. NCG: diagonal preconditioner
    yik_sos = sum(conj(smap) .* ydata; dims=3)[:,:,1,:]; # coil combine => (nY, nX, 1, nCon) => (nY, nX, nCon)

    (yik_sos_scaled, scale) = b0scale(yik_sos, TE); # fig = plt_images(permutedims(abs.(yik_sos_scaled), [3,1,2]),width=10, height=5)

    # b0 init
    # finit = b0init(ydata, TE; smap); plt_image(finit; title="b0init")

    yik_scale = ydata / scale;
    fmap_run = (niter, precon, track; kwargs...) -> b0map(yik_scale, TE; smap, mask, order=1, l2b=l2b, gamma_type=:PR, niter, precon, track, kwargs...)

    function runner(niter, precon; kwargs...)
        (fmap, times, out) = fmap_run(niter, precon, true; kwargs...) # tracking run
        return (fmap, out.fhats, out.costs, times)
    end;

    (fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter, precon);
    b0 = fmap_cg_d;  # Hz
    return b0, yik_sos
end

function get_mask(image; threshold=0)
    image = abs.(image);
    mask = zeros(size(image));
    mask[image.>threshold] .= 1;
    mask = isone.(mask);
    return mask
end
