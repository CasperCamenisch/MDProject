# Molecular Dynamics Simulation

# Loading packages
using Distributions
using Plots
using StatsPlots
using LinearAlgebra

# Loading code files
include("thermostats.jl")
include("particles.jl")
include("settings_reader.jl")
include("input_reader.jl")
include("input_generator.jl")
include("output_writer.jl")

# Periodic distance calculations
function dist_vector_periodic(a_pos::Vector{Float32}, b_pos::Vector{Float32}, L::Float32, dim::Int8)
    vect::Vector{Float32} = [0.0 for _ in 1:dim]
    for n in 1:dim
        dx::Float32 = a_pos[n] - b_pos[n]
        if abs(dx) < 0.5*L
            vect[n] = dx
        else
            if dx > 0
                vect[n] = -(L-dx)
            else
                vect[n] = L+dx
            end
        end
    end
    return vect
end

function sq_dist_periodic(a_pos::Vector{Float32}, b_pos::Vector{Float32}, L::Float32, dim::Int8)
    s::Float32 = 0.0

    for n in 1:dim
        dx::Float32 = abs(a_pos[n] - b_pos[n])
        if dx > 0.5*L
            s += (L - dx)^2
        else
            s += dx^2
        end
    end

    return s
end

function boundary_correction(x::Float32, BoxSize::Float32)
    if x < 0.0
        x = BoxSize + (x % BoxSize)
    elseif x > BoxSize
        x = x % BoxSize
    end

    return x
end

function correct_positions(particles::Particles, BoxSize::Float32)
    # Corrects postitions if over the Box boundary
    particles.position = [[boundary_correction(x, BoxSize) for x in pos] for pos in particles.position]
end

# Verlet neighbour list
function neighbour_list(particles::Particles, N::Int32, BoxSize::Float32, dim::Int8, cutoff::Float32, r_l::Float32)

    r_list_cutoff::Float32 = cutoff + r_l
    clusters::Vector{Vector{Int32}} = []       # Stores lists of neighbours of all particles

    for i in 1:N
        n_list::Vector{Int32} = []     # Stores indices of neighbours of given particle i
        for j in i+1:N

            r_sq::Float32 = sq_dist_periodic(particles.position[i], particles.position[j], BoxSize, dim)

            if r_sq < r_list_cutoff^2
                push!(n_list, j)
            end
        end
        push!(clusters, n_list)
    end

    return clusters
end

# VdW force calculation
function calculate_forces(N::Int32, particles::Particles, cutoff::Float32, BoxSize::Float32, dim::Int8, clusters::Vector{Vector{Int32}})
    accelerations::Vector{Vector{Float32}} = [[0.0 for _ in 1:dim] for _ in 1:N]
    virial::Float32 = 0.0
    n_start::Int32 = 1
    # Loop over all clusters
    for i in 1:N-1
        for j in clusters[i]
            r_sq::Float32 = sq_dist_periodic(particles.position[i], particles.position[j], BoxSize, dim)
            if r_sq < cutoff^2
                # Get LJ parameter values for interaction
                sigma::Float32 = LJ_sigma_parameters[particles.type[i]][particles.type[j]]
                epsilon::Float32 = LJ_epsilon_parameters[particles.type[i]][particles.type[j]]
                # Intermediate calculation steps
                sigr_2 = (sigma^2)/r_sq
                sigr_6 = sigr_2^3
                sigr_12 = sigr_6^2
                # Calculate potential gradient
                d_potential::Float32 = epsilon * 24.0 * sigr_2 * (2.0 * sigr_12 - sigr_6)
                r::Float32 = sqrt(r_sq)
                # Virial coefficient
                virial += d_potential * r
                # Update accelerations
                dist::Vector{Float32} = normalize(dist_vector_periodic(particles.position[i], particles.position[j], BoxSize, dim))
                a::Vector{Float32} = d_potential * dist / sqrt(r_sq)
                accelerations[i] += a
                accelerations[j] -= a
            end
        end
    end

    # Adjust accelerations to particle masses
    for i in 1:N
        accelerations[i] = accelerations[i] / particles.mass[i]
    end

    return accelerations, virial
end

# Time integration
function velocity_verlet(h::Float32, N::Int32, particles::Particles, target_temperature::Float32, thermo_type::String, cutoff::Float32, BoxSize::Float32, dim::Int8, a_old::Vector{Vector{Float32}}, clusters::Vector{Vector{Int32}})
    # Returns coordinates and velocities after 1 time-step
    particles.position += particles.velocity*h + 0.5*a_old*h^2
    a_new::Vector{Vector{Float32}}, virial::Float32 = calculate_forces(N, particles, cutoff, BoxSize, dim, clusters)
    particles.velocity += 0.5*(a_old + a_new)*h
    # Thermostat
    thermostat(thermo_type, N, particles, BoxSize, dim, target_temperature, h)

    return a_new
end


# Solve ODE
function odeSolver(h, nSteps, N, particles, target_temperature, thermo_type, cutoff, BoxSize, dim, interval, density, V, r_l, outfile)
    
    # Calculate neighbour list update steps
    v::Vector{Float32} = [0.0 for _ in 1:N]
    for i in 1:N
        v[i] = sqrt(dot(particles.velocity[i], particles.velocity[i]))
    end
    v_max::Float32 = maximum(v)
    list_interval::Int32 = floor(0.5 * r_l / (v_max*h))
    clusters::Vector{Vector{Int32}} = neighbour_list(particles, N, BoxSize, dim, cutoff, r_l)
    
    # Get inital acceleration
    accel::Vector{Vector{Float32}}, virial::Float32 = calculate_forces(N, particles, cutoff, BoxSize, dim, clusters)
    
    # Simulation loop
    for i in 1:nSteps
        accel = velocity_verlet(h, N, particles, target_temperature, thermo_type, cutoff, BoxSize, dim, accel, clusters)
        correct_positions(particles, BoxSize)
        e_kin::Float32, temperature::Float32 = calculate_temperature(N, particles, BoxSize, dim)
        pressure::Float32 = density * temperature + virial/V

        if i % interval == 0
            write_step(outfile, i, N, particles, temperature, pressure)
        end

        if i % list_interval == 0
            
            # Update neighbour list
            clusters = neighbour_list(particles, N, BoxSize, dim, cutoff, r_l)
            
            # Set new interval
            for i in 1:N
                v[i] = sqrt(dot(particles.velocity[i], particles.velocity[i]))
            end
            v_max = maximum(v)
            list_interval = floor(0.5 * r_l / (v_max*h))
            if list_interval == 0
                list_interval = 1
            end
        end
    end

    return 0
end

# Main fucntion
function main()
    println("Initializing...")
    settings::SETTINGS = read_settings(ARGS[1])
    # Testing box parameters
    N::Int32 = settings.N
    BoxSize::Float32 = settings.L
    dim::Int8 = settings.DIM
    V::Float32 = BoxSize^dim
    density::Float32 = N/V
    # Particles
    if settings.INPUT_FILE != ""
        particles = read_data(settings.INPUT_FILE)
    else
        particles = generate_grid(N, settings.TYPE, BoxSize, dim)
    end
    target_temperature::Float32 = settings.T
    thermo_type::String = settings.THERMO

    # Set initial particle velocities
    temp_init(N, particles, dim, target_temperature)

    # Calculate initial temperature
    init_e::Float32, init_temp::Float32 = calculate_temperature(N, particles, BoxSize, dim)

    # Solver settings
    cutoff::Float32 = settings.CUTOFF
    h::Float32 = settings.STEPSIZE
    steps::Int32 = settings.STEPS
    r_l::Float32 = settings.LIST_RADIUS
    interval::Int32 = settings.INTERVAL

    # Running simulation
    println("Solver started")
    write_start(settings.OUTPUT_FILE, h, N, particles, thermo_type, target_temperature, init_temp)
    @time odeSolver(h, steps, N, particles, target_temperature, thermo_type, cutoff, BoxSize, dim, interval, density, V, r_l, settings.OUTPUT_FILE)
    println("Done!")
end

main()