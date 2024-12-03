export normalOperator

function LinearOperatorCollection.NormalOpi1(::Type{T};
  parent, weights, tmp::Union{Nothing,Vector{T}}) where T <: Number
  if tmp == nothing
    return NormalOpi1Impl(parent, weights)
  else
    return NormalOpi1Impl(parent, weights, tmp)
  end
end

mutable struct NormalOpi1Impl{T,S,D,V} <: NormalOpi1{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod! :: Function
  tprod! :: Nothing
  ctprod! :: Nothing
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  args5 :: Bool
  use_prod5! :: Bool
  allocated5 :: Bool
  Mv5 :: Vector{T}
  Mtu5 :: Vector{T}
  parent::S
  weights::D
  tmp::V
end

LinearOperators.storage_type(op::NormalOpi1Impl) = typeof(op.Mv5)

function NormalOpi1Impl(parent, weights)
  T = promote_type(eltype(parent), eltype(weights))
  tmp = Vector{T}(undef, size(parent, 1))
  return NormalOpi1Impl(parent, weights, tmp)
end

function NormalOpi1Impl(parent, weights, tmp::Vector{T}) where T

  function produ!(y, parent, tmp, x)
    mul!(tmp, parent, x)
    mul!(tmp, weights, tmp) # This can be dangerous. We might need to create two tmp vectors
    return mul!(y, adjoint(parent), tmp)
  end

  return NormalOpi1Impl(size(parent,2), size(parent,2), false, false
         , (res,x) -> produ!(res, parent, tmp, x)
         , nothing
         , nothing
         , 0, 0, 0, false, false, false, T[], T[]
         , parent, weights, tmp)
end

function Base.copy(S::NormalOpi1Impl)
  return NormalOpi1Impl(copy(S.parent), S.weights, copy(S.tmp))
end

function normalOperator(parent, weights=opEye(eltype(parent), size(parent,1)))
  return NormalOpi1Impl(parent, weights)
end

