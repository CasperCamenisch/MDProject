# Generating input particle read_data
include("particles.jl")

function generate_data(N, type, BoxSize, dim, temperature)
    atom_types = [type for _ in 1:N]
    if type == "Ar"
        masses = [39.962 for _ in 1:N]
        charges = [0 for _ in 1:N]
    end
    r_array = [rand(Uniform(0.0, BoxSize), dim) for _ in 1:N]
    v_array = [rand(Uniform(-0.5, 0.5), dim) for _ in 1:N]  # TODO: adjust to reflect temperature
    v_array = [[0.0 for _ in 1:dim] for _ in 1:N]

    particles = Particles(atom_types, masses, charges, r_array, v_array)

    return particles
end

function generate_grid(N, type, BoxSize, dim)
    # ONLY works if N^(1/dim) is a whole number -> 1, 8, 27, 64, 125, 216, 343, 512, 729, 1000, 2197, 4096, 8000, 9261, 10648, 15625 
    # used as such: p = generate_grid(28, "Ar", 3000, 3)
    atom_types = [type for _ in 1:N]
    if type == "Ar"
        masses = [39.962 for _ in 1:N]
        charges = [0 for _ in 1:N]
    end
    if type == "Kr"
        masses = [83.80 for _ in 1:N]
        charges = [0 for _ in 1:N]
    end
    if type == "Xe"
        masses = [131.30 for _ in 1:N]
        charges = [0 for _ in 1:N]
    end

    n_row::Int = round(N^(1/dim))
    stepsize = BoxSize/n_row

    if dim == 2
        r_array = [[i*stepsize, j*stepsize] for i in 0:n_row-1 for j in 0:n_row-1]
        println("N = ", length(r_array))
        println("interparticel distance = $stepsize")
    elseif dim == 3
        r_array = [[i*stepsize, j*stepsize, k*stepsize] for i in 0:n_row-1 for j in 0:n_row-1 for k in 0:n_row-1]
    end

    v_array = [[0.0 for _ in 1:dim] for _ in 1:N]

    particles = Particles(atom_types, masses, charges, r_array, v_array)

    return particles
end

function N2_system(type, BoxSize)
    atom_types = [type for _ in 1:2]
    if type == "Ar"
        masses = [39.962 for _ in 1:2]
        charges = [0 for _ in 1:2]
    end
    r_array = [[0.5*BoxSize,0.5*BoxSize - 300.0], [0.5*BoxSize,0.5*BoxSize + 300.0]]
    v_array = [[0.0,0.0], [0.0,0.0]]

    particles = Particles(atom_types, masses, charges, r_array, v_array)

    return particles
end