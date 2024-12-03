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
    x = solve(solver, kdata)
    x = reshape(x, recParams[:reconSize])
    return x
end