function recon_HOOp(HOOp::HOOp{Complex{T}}, Data::AbstractArray{Complex{T},2}, recParams::Dict) where T<:AbstractFloat
    recoParams = merge(defaultRecoParams(), recParams)
    reg = Regularization(recParams[:regularization], recParams[:λ]; shape=recParams[:reconSize])

    nSample, nCha = size(Data)
    Data = vec(Data)
    E = HOOp
    EᴴE = normalOperator(E)
    solver = createLinearSolver(recParams[:solver], E; AᴴA=EᴴE, reg=reg, recoParams...)
    x = solve(solver, Data; recoParams...)
    x = reshape(x, recParams[:reconSize])
    return x
end

function recon_HOOp(HOOp::HOOp{Complex{T}}, Data::AbstractArray{Complex{T},2}, weight::AbstractVector{Complex{T}}, recParams::Dict) where T<:AbstractFloat
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
