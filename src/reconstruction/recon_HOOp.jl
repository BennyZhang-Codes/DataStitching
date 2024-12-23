function recon_HOOp(HOOp::HighOrderOp, acqData::AcquisitionData, recParams::Dict) :: Matrix
    recoParams = merge(defaultRecoParams(), recParams)
    numContr, numChan = MRIReco.numContrasts(acqData), MRIReco.numChannels(acqData);
    reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, recParams);
    senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan);
    smaps = senseMaps[:,:,1,:];

    kdata = multiCoilData(acqData, 1, 1, rep=1) .* repeat(weights[1], numChan)
    W = WeightingOp(ComplexF64; weights=weights[1], rep=numChan)
    E = ∘(W, HOOp, isWeighting=false)
    EᴴE = normalOperator(E)
    solver = createLinearSolver(recParams[:solver], E; AᴴA=EᴴE, recoParams...)
    x = solve(solver, kdata)
    x = reshape(x, recParams[:reconSize])
    return x
end


function recon_HOOp(HOOp::HighOrderOp_i2, acqData::AcquisitionData, recParams::Dict) :: Matrix
    recoParams = merge(defaultRecoParams(), recParams)
    numContr, numChan = MRIReco.numContrasts(acqData), MRIReco.numChannels(acqData);
    reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, recParams);
    senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan);
    smaps = senseMaps[:,:,1,:];

    kdata = multiCoilData(acqData, 1, 1, rep=1) .* repeat(weights[1], numChan)
    W = WeightingOp(ComplexF64; weights=weights[1], rep=numChan)
    E = ∘(W, HOOp, isWeighting=false)
    EᴴE = normalOperator(E)
    solver = createLinearSolver(recParams[:solver], E; AᴴA=EᴴE, reg=reg, recoParams...)
    x = solve(solver, kdata; recoParams...)
    x = reshape(x, recParams[:reconSize])
    return x
end


function recon_HOOp(HOOp::HighOrderOpv2{Complex{T}}, Data::AbstractArray{Complex{T},2}, weight::AbstractVector{Complex{T}}, recParams::Dict) where T<:AbstractFloat
    recoParams = merge(defaultRecoParams(), recParams)
    reg = Regularization(recParams[:regularization], recParams[:λ]; shape=recParams[:reconSize])

    nSample, nCha = size(Data)
    Data = vec(Data) .* repeat(weight, nCha)
    W = WeightingOp(Complex{T}; weights=weight, rep=nCha)
    E = ∘(W, HOOp, isWeighting=false)
    EᴴE = normalOperator(E)
    solver = createLinearSolver(recParams[:solver], E; AᴴA=EᴴE, reg=reg, recoParams...)
    x = solve(solver, Data; recoParams...)
    x = reshape(x, recParams[:reconSize])
    return x
end

function recon_HOOp(HOOp::HighOrderOpv2_i2{Complex{T}}, Data::AbstractArray{Complex{T},2}, weight::AbstractVector{Complex{T}}, recParams::Dict) where T<:AbstractFloat
    recoParams = merge(defaultRecoParams(), recParams)
    reg = Regularization(recParams[:regularization], recParams[:λ]; shape=recParams[:reconSize])

    nSample, nCha = size(Data)
    Data = vec(Data) .* repeat(weight, nCha)
    W = WeightingOp(Complex{T}; weights=weight, rep=nCha)
    E = ∘(W, HOOp, isWeighting=false)
    EᴴE = normalOperator(E)
    solver = createLinearSolver(recParams[:solver], E; AᴴA=EᴴE, reg=reg, recoParams...)
    x = solve(solver, Data; recoParams...)
    x = reshape(x, recParams[:reconSize])
    return x
end