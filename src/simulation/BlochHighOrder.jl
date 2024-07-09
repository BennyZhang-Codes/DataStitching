
Base.@kwdef struct BlochHighOrder <: SimulationMethod 
    ho0::Bool = true
    ho1::Bool = true
    ho2::Bool = true
    Δw_excitation::Bool = true
    Δw_precession::Bool = true
    name::String = string(Int(ho0)) * string(Int(ho1)) * string(Int(ho2))
end

function BlochHighOrder(name::String)
    @assert name == "000" || name == "100" || name == "010" || name == "001" || name == "110" || name == "101" || 
    name == "011" || name == "111" "Invalid name for BlochHighOrder simulation method. Valid names are: 000, 100, 010, 001, 110, 101, 011, 111"
    ho0 = Bool(parse(Int64, name[1]))
    ho1 = Bool(parse(Int64, name[2]))
    ho2 = Bool(parse(Int64, name[3]))
    return BlochHighOrder(ho0=ho0, ho1=ho1, ho2=ho2)
end

function BlochHighOrder(name::String, Δw_excitation::Bool, Δw_precession::Bool)
    @assert name == "000" || name == "100" || name == "010" || name == "001" || name == "110" || name == "101" || 
    name == "011" || name == "111" "Invalid name for BlochHighOrder simulation method. Valid names are: 000, 100, 010, 001, 110, 101, 011, 111"
    ho0 = Bool(parse(Int64, name[1]))
    ho1 = Bool(parse(Int64, name[2]))
    ho2 = Bool(parse(Int64, name[3]))
    return BlochHighOrder(ho0=ho0, ho1=ho1, ho2=ho2, Δw_excitation=Δw_excitation, Δw_precession=Δw_precession)
end

Base.show(io::IO, b::BlochHighOrder) = begin
	print(io, "BlochHighOrder[[$(b.name)] zero order=$(b.ho0) | first order=$(b.ho1) | second order=$(b.ho2) | Δw_excitation=$(b.Δw_excitation) | Δw_precession=$(b.Δw_precession)]")
end