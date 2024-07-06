struct sphericalharmonic
    name::String
    unit::String
    cumtrapz_unit::String
    expression::String
end
Base.show(io::IO, s::sphericalharmonic) = print(io, s.name, " [", s.unit, "] = ", s.expression)

Base.@kwdef struct SphericalHarmonics
    h0::sphericalharmonic =sphericalharmonic("h0",    "T",    "", "1")
    h1::sphericalharmonic =sphericalharmonic("h1",  "T/m", "m⁻¹", "x")
    h2::sphericalharmonic =sphericalharmonic("h2",  "T/m", "m⁻¹", "y")
    h3::sphericalharmonic =sphericalharmonic("h3",  "T/m", "m⁻¹", "z")
    h4::sphericalharmonic =sphericalharmonic("h4", "T/m²", "m⁻²", "xy")
    h5::sphericalharmonic =sphericalharmonic("h5", "T/m²", "m⁻²", "zy")
    h6::sphericalharmonic =sphericalharmonic("h6", "T/m²", "m⁻²", "3z² - (x² + y² + z²)")
    h7::sphericalharmonic =sphericalharmonic("h7", "T/m²", "m⁻²", "xz")
    h8::sphericalharmonic =sphericalharmonic("h8", "T/m²", "m⁻²", "x² - y²")
    dict::Dict =Dict("h0"=>h0, "h1"=>h1, "h2"=>h2, "h3"=>h3, "h4"=>h4, "h5"=>h5, "h6"=>h6, "h7"=>h7, "h8"=>h8)
end

Base.show(io::IO, s::SphericalHarmonics) = begin
    print(io, "Spherical Harmonics")
    print(io, "\n  ", s.h0)
    print(io, "\n  ", s.h1)
    print(io, "\n  ", s.h2)
    print(io, "\n  ", s.h3)
    print(io, "\n  ", s.h4)
    print(io, "\n  ", s.h5)
    print(io, "\n  ", s.h6)
    print(io, "\n  ", s.h7)
    print(io, "\n  ", s.h8)
end