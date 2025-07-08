"""
    BHO = BlochHighOrder(ho0, ho1, ho2, ho3, Δw_excitation, Δw_precession)

a Simulation Method, which defines some parameters for the simulation with dynamic fields (up to second order).

# Arguments
- `ho0`: (`::Bool`), zeroth order
- `ho1`: (`::Bool`), first order
- `ho2`: (`::Bool`), second order
- `ho3`: (`::Bool`), third order
- `Δw_excitation`: (`::Bool`), whether to include ΔB₀ in excitation
- `Δw_precession`: (`::Bool`), whether to include ΔB₀ in precession

# Returns
- `BHO`: (`::BlochHighOrder`) BlochHighOrder struct
"""
Base.@kwdef struct BlochHighOrder <: SimulationMethod 
    ho0::Bool = true
    ho1::Bool = true
    ho2::Bool = true
    ho3::Bool = true
    Δw_excitation::Bool = true
    Δw_precession::Bool = true
    name::String = string(Int(ho0)) * string(Int(ho1)) * string(Int(ho2)) * string(Int(ho3))
end

function BlochHighOrder(name::String)
    @assert name == "0000" || name == "1000" || name == "0100" || name == "0010" || name == "0001" ||
    name == "1100" || name == "1001" || name == "0011" || name == "1010" || name == "0110" || name == "0101" ||
    name == "1110" || name == "1101" || name == "1011" || name == "0111" ||
    name == "1111" "Invalid name for BlochHighOrder simulation method. Valid names are: 0000, 1000, 0100, 0010, 0001, 1100, 1001, 0011, 1010, 0110, 0101, 1110, 1101, 1011, 0111, 1111"
    ho0 = Bool(parse(Int64, name[1]))
    ho1 = Bool(parse(Int64, name[2]))
    ho2 = Bool(parse(Int64, name[3]))
    ho3 = Bool(parse(Int64, name[4]))
    return BlochHighOrder(ho0=ho0, ho1=ho1, ho2=ho2, ho3=ho3)
end

function BlochHighOrder(name::String, Δw_excitation::Bool, Δw_precession::Bool)
    @assert name == "0000" || name == "1000" || name == "0100" || name == "0010" || name == "0001" ||
    name == "1100" || name == "1001" || name == "0011" || name == "1010" || name == "0110" || name == "0101" ||
    name == "1110" || name == "1101" || name == "1011" || name == "0111" ||
    name == "1111" "Invalid name for BlochHighOrder simulation method. Valid names are: 0000, 1000, 0100, 0010, 0001, 1100, 1001, 0011, 1010, 0110, 0101, 1110, 1101, 1011, 0111, 1111"
    ho0 = Bool(parse(Int64, name[1]))
    ho1 = Bool(parse(Int64, name[2]))
    ho2 = Bool(parse(Int64, name[3]))
    ho3 = Bool(parse(Int64, name[4]))
    return BlochHighOrder(ho0=ho0, ho1=ho1, ho2=ho2, ho3=ho3, Δw_excitation=Δw_excitation, Δw_precession=Δw_precession)
end

Base.show(io::IO, b::BlochHighOrder) = begin
	print(io, "BlochHighOrder[[$(b.name)] zero order=$(b.ho0) | first order=$(b.ho1) | second order=$(b.ho2) | third order=$(b.ho3) | Δw_excitation=$(b.Δw_excitation) | Δw_precession=$(b.Δw_precession)]")
end