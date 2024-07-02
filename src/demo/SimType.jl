
Base.@kwdef struct SimType
    B0::Bool = true
    T2::Bool = true
    ss::Int64 = 3
    name::String = "B0$(B0 ? "w" : "wo")_T2$(T2 ? "w" : "wo")_ss$(ss)"
end

Base.show(io::IO, b::SimType) = begin
	print(io, "SimType[[$(b.name)] B0=$(b.B0) | T2=$(b.T2) | ss=$(b.ss) ]")
end

function SimType(name::String)
    @assert length(split(name,"_")) == 3 "Invalid name for SimType. Valid name like: B0wo_T2w_ss3"
    B0, T2, ss = split(name,"_")
    @assert !isnothing(findfirst("B0", B0)) && !isnothing(findfirst("T2", T2)) && !isnothing(findfirst("ss", ss)) "Invalid name for SimType. Valid name like: B0wo_T2w_ss3"
    B0, T2, ss = split(B0, "B0")[2], split(T2, "T2")[2], split(ss, "ss")[2]
    B0, T2, ss = B0 == "w", T2 == "w", parse(Int, ss)
    return SimType(B0=B0, T2=T2, ss=ss)
end
