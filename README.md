# Instructions

The usage of this program requires a working JULIA environment.

## Settings files
Settings files are text files that provide the setting of an MD simulation. They have to be given as an argument as in
```
julia MD_sim.jl my_settings.txt
```
There are two examples of settings files given in the /src/ directory.
- "settings\_example1.txt" is an example of a settings file that where the initial configuration is generated to be a regular grid of a given atom.
- "setting\_example2.txt" is an example fo a settings file that takes the initial configuration from an XYZ file that has to be provided.

## Parameters

- N | the number of particles.
- L | the sidelength of the simulation box (in pm).
- DIM | the number of spacial dimensions. Usually 3 but 2 can be used for testing or to get nice visualizations.
- INPUT\_FILE | The file providing the initial configuration of the system in XYZ format. Not necessary if you want to regular grid as the initial configuration.
- TYPE | Type of the atoms. Not necessary if an input file is given. Options are "Ar", "Kr" or "Xe".
- T | the initial temperature and also the target temperature (in K) if a thermostat is used.
- THERMO | The thermostat. Possible options are "None", "Rescale", "Berendsen" or "Andersen".
- CUTOFF | The cutoff radius (in pm) for the LJ potential.
- LIST\_RADIUS | The *additional* radius for creating the Verlet neighbour list.
- STEPSIZE | Length of a single time-step in the simulation (in fs).
- STEPS | Number of time-steps for the simulation.
- INTERVAL | Number of time-steps that are skipped when the output is written to a file.
- OUTPUT\_FILE | Name of the file where the output is written.

