# Read input file (XYZ format)
include("particles.jl")

function read_data(filename)
    atom_types = []
    masses = []
    charges = []
    r_array = []
    v_array = []

    
    f = open(filename) do f
        readline(f) # Skipping the first 2 header lines
        readline(f)

        while ! eof(f)
            line = split(readline(f))
            type = line[1]
            x = parse(Float32, line[2]) * 100.0
            y = parse(Float32, line[3]) * 100.0
            z = parse(Float32, line[4]) * 100.0
            push!(atom_types, type)
            push!(r_array, [x, y, z])
            push!(v_array, [0.0, 0.0, 0.0])
            push!(charges, 0)
            
            if type == "Ar"
                push!(masses, 39.962)
            end
        end
    end

    # Min / Max values for each dimension
    #println(reduce((x,y) -> max.(x,y), r_array))
    #println(reduce((x,y) -> min.(x,y), r_array))

    particles = Particles(atom_types, masses, charges, r_array, v_array)

    return particles #[atom_types, x_array, y_array, z_array]
end

#= 
Used as such:
p = read_data("argon.xyz")
=#