# Class and Function definitions of thermostats for use in MD Simulation
include("particles.jl")

function calculate_temperature(N::Int32, particles::Particles, BoxSize::Float32, dim::Int8)
    e_kin::Float32 = 0.0

    for i in 1:N
        velocity_sq::Float32 = dot(particles.velocity[i], particles.velocity[i])
        e_kin += 0.5 * particles.mass[i] * velocity_sq
    end

    e_kin_avg::Float32 = e_kin/N
    reverse_kb::Float32 = 120.33173051468289 #0.724297051603992  # 1/boltzmann constant
    temperature::Float32 = 2.0/dim * e_kin_avg * reverse_kb

    return e_kin_avg, temperature
end

function velocity_rescale(N::Int32, particles::Particles, BoxSize::Float32, dim::Int8, target_temperature::Float32)
    """ 
    Rescales velocities ensuring that T is constant
    """
    _, temp::Float32 = calculate_temperature(N, particles, BoxSize, dim)
    scaling::Float32 = sqrt(target_temperature/temp)
    if scaling == 0
        println(temp)
    end
    # Update velocities
    particles.velocity = scaling * particles.velocity
end

function Berendsen(N::Int32, particles::Particles, BoxSize::Float32, dim::Int8, target_temperature::Float32, h::Float32)
    """ 
    Rescales velocitities ensuring that T is moved towards the target_temperature
    """
    tau::Float32 = 100.0 * h     # if tau = 1*h we have the simple velocity rescale method
    _, temp::Float32 = calculate_temperature(N, particles, BoxSize, dim)
    scaling::Float32 = sqrt(1 + h/tau * (target_temperature/temp - 1))       # due to number of dimension the scaling might be off
    # Update velocities
    particles.velocity = scaling * particles.velocity
end

function Andersen(N::Int32, particles::Particles, dim::Int8, target_temperature::Float32, dt::Float32)
    """ 
    Particle collisions at given frequency with heat bath to maintain constant T
    """
    nu::Float32 = 0.1 #0.1        # select 10% of particles
    conv_dim::Float32 = 1 #1/sqrt(3)  #(1/dim)^(1/2.72) #dim^(1/3)  #dim^(-0.5)  # convert to velocity to correct dimensions
    kb::Float32 = 0.00831036 #1.380649   # Boltzmann constant
    for i in 1:N
        random_draw::Float32 = rand(Uniform(0.0, 1.0)) # draw from uniform distribution
        if random_draw < nu*dt  # dt = h
            # prefactor
            a::Float32 = conv_dim*sqrt(kb * target_temperature / particles.mass[i])
            # draw from Maxwell-Boltzmann distribution at target_temperature
            particles.velocity[i] = [rand([1, -1])*a*rand(Chi(1)) for _ in 1:dim]
        end
    end
end

function temp_init(N::Int32, particles::Particles, dim::Int8, target_temperature::Float32)
    """
    Sets particle velocities, taken from Maxwell-Boltzmann distribution 
    """
    # Set initial particle velocities
    conv_dim::Float32 = 1 #1/(sqrt(dim))    # convert to velocity to correct dimensions
    kb::Float32 = 0.00831036 #1.380649   # Boltzmann constant
    for i in 1:N
        # prefactor from
        a::Float32 = conv_dim*sqrt(kb * target_temperature / particles.mass[i])
        # draw from Maxwell-Boltzmann distribution at target_temperature
        particles.velocity[i] = [rand([1, -1])*a*rand(Chi(1)) for _ in 1:dim]
    end
end

function thermostat(type::String, N::Int32, particles::Particles, BoxSize::Float32, dim::Int8, target_temperature::Float32, h::Float32)
    if type == "Rescale"
        temp = velocity_rescale(N, particles, BoxSize, dim, target_temperature)
    elseif type == "Berendsen"
        temp = Berendsen(N, particles, BoxSize, dim, target_temperature, h)
    elseif type == "Andersen"
        temp = Andersen(N, particles, dim, target_temperature, h)
    elseif type == "None"
        # Does nothing
    else
        println("Thermostat not defined!")
        exit()
    end
end