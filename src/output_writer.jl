# Write the output to an .out file

function write_start(filename, stepsize, N, particles, thermostat, target_temperature, T)
    file = open(filename, "w")
    write(file, "\nStepsize =  $stepsize fs\nThermostat = $thermostat\nTarget Temperature = $target_temperature K\nStarting Temperature = $T K\n")
    for i in 1:N
        type = particles.type[i]
        x = particles.position[i][1] * 0.01
        y = particles.position[i][2] * 0.01
        z = particles.position[i][3] * 0.01
        write(file, "$type\t$x\t$y\t$z\n")
    end
    close(file)
end

function write_step(filename, i_step, N, particles, T, pressure)
    file = open(filename, "a")
    write(file, "\nstep $i_step\nT = $T\nP = $pressure\n")
    for i in 1:N
        type = particles.type[i]
        x = particles.position[i][1] * 0.01
        y = particles.position[i][2] * 0.01
        z = particles.position[i][3] * 0.01
        write(file, "$type\t$x\t$y\t$z\n")
    end
    close(file)

end