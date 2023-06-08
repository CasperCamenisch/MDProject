# Struct for particles (type, coords, charge etc.)

mutable struct Particles
    type::Vector{String}
    mass::Vector{Float32}
    charge::Vector{Float32}
    position::Vector{Vector{Float32}}
    velocity::Vector{Vector{Float32}}
end


# Lookup table for parameters (epsilon, sigma)
# Data from Liquid-State Physical Chemistry, First Edition. Gijsbertus de With.
#© 2013 Wiley-VCH Verlag GmbH & Co. KGaA. Published 2013 by Wiley-VCH Verlag GmbH & Co. KGaA.

# Epsilon in epsilon/kb (K) units converted to pm^2 amu fs^-2 by the boltzmann constant
#kb = 8.31036
kb::Float32 = 0.00831036
LJ_epsilon_parameters = Dict(
"Ar" => Dict(
    "Ar" => kb * 111.84,
    "Kr" => kb * sqrt(111.84 * 154.87),
    "Xe" => kb * sqrt(111.84 * 213.96)
    ),
"Kr" => Dict(
    "Ar" => kb * sqrt(111.84 * 154.87),
    "Kr" => kb * 154.87,
    "Xe" => kb * sqrt(154.87 * 213.96)
    ),
"Xe" => Dict(
    "Ar" => kb * sqrt(111.84 * 213.96),
    "Kr" => kb * sqrt(154.87 * 213.96),
    "Xe" => kb * 213.96
    )
)

# Sigma in pm (10^-12 m)
LJ_sigma_parameters = Dict(
"Ar" => Dict(
    "Ar" => 362.3,
    "Kr" => 0.5*(362.3 + 389.5),
    "Xe" => 0.5*(362.3 + 426.0)
    ),
"Kr" => Dict(
    "Ar" => 0.5*(362.3 + 389.5),
    "Kr" => 389.5,
    "Xe" => 0.5*(389.5 + 426.0)
    ),
"Xe" => Dict(
    "Ar" => 0.5*(362.3 + 426.0),
    "Kr" => 0.5*(389.5 + 426.0),
    "Xe" => 426.0
    )
)