
function get_mask(image; threshold=0)
    image = abs.(image);
    mask = zeros(size(image));
    mask[image.>threshold] .= 1;
    mask = isone.(mask);
    return mask
end

function read_gre(seqfile::String, mrdfile::String) :: Tuple
    @assert isfile(seqfile) "sequence file not found: $seqfile"
    @assert isfile(mrdfile) "raw data file not found: $mrdfile"

    # read sequence and trajectory
    seq = read_seq(seqfile);
    _, ktraj = get_kspace(seq);
    TE  = seq.DEF["TE"];  # s
    ReadoutMode = seq.DEF["ReadoutMode"]; # "Bipolar" or "Monopolar"

    # get kspace data and shape
    raw  = RawAcquisitionData(ISMRMRDFile(mrdfile))
    shape = get_ksize(raw);
    nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape;
    kdata = get_kdata(raw, shape);
    kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
    kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];
    @info "read multi-echo gre data\n  size=$(size(kdata))\n  dims=$(kdims)\n  ReadoutMode=$(ReadoutMode)"
    if ReadoutMode=="Bipolar"
        kdata[:,:,:, collect(range(2, nCon, step=2))] = kdata[:,:,end:-1:1, collect(range(2, nCon, step=2))]; # reverse the even echoes because of bipolar readout
    end

    # prepare acqData for coil sensitivity estimation
    Con = 1
    ktr = reshape(ktraj[:,1:2], nX, nCon*nY, :)[:, collect(range(Con,nCon*nY,step=nCon)), :];
    ktr = reshape(ktr, nX*nY, :);
    tr  = Trajectory(T.(ktr'), nY, nX, circular=false, cartesian=true); # plot_traj2d(tr)

    k_img = Array{Complex{T},2}(undef,nX*nY,nCha);
    for y = 1:nY
        k_img[(y-1)*nX+1:y*nX,:] = transpose(kdata[:, y, : ,Con])
    end
    # plt_image(abs.(kdata[1, :, : ,Con]).^0.02)
    # fig, ax = plt.subplots(1,1); ax.plot(abs.(k_img[:,1]));
    dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
    dat[1,1,1]   = reshape(k_img,:,nCha);
    acqData      = AcquisitionData(tr, dat, encodingSize=(nX, nY));

    return kdata, ktraj, kdims, shape, TE, ReadoutMode, acqData, raw
end

using MRIFieldmaps: b0scale, b0init, b0map
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
    b0 = fmap_cg_d * -1;  # Hz
    return b0, yik_sos
end