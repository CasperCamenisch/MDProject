# Parse simulation settings

struct SETTINGS
    N::Int32
    DIM::Int8
    L::Float32
    INPUT_FILE::String
    TYPE::String
    THERMO::String
    T::Float32
    CUTOFF::Float32
    STEPSIZE::Float32
    STEPS::Int32
    LIST_RADIUS::Float32
    INTERVAL::Int32
    OUTPUT_FILE::String
end

function read_settings(filename)
    N::Int32 = 0
    DIM::Int8 = 0
    L::Float32 = 0
    INPUT_FILE::String = ""
    TYPE::String = ""
    THERMO::String = ""
    T::Float32 = 0
    CUTOFF::Float32 = 0
    STEPSIZE::Float32 = 0
    STEPS::Int32 = 0
    LIST_RADIUS::Float32 = 0
    INTERVAL::Int32 = 0
    OUTPUT_FILE::String = ""

    f = open(filename) do f

        while ! eof(f)
            line = split(readline(f))
            parameter = line[1]
            setting = line[2]

            if parameter == "N"
                N = parse(Int32, setting)
            end
            if parameter == "DIM"
                DIM = parse(Int8, setting)
            end
            if parameter == "L"
                L = parse(Float32, setting)
            end
            if parameter == "INPUT_FILE"
                INPUT_FILE = setting
            end
            if parameter == "TYPE"
                TYPE = setting
            end
            if parameter == "THERMO"
                THERMO = setting
            end
            if parameter == "T"
                T = parse(Float32, setting)
            end
            if parameter == "CUTOFF"
                CUTOFF = parse(Float32, setting)
            end
            if parameter == "STEPSIZE"
                STEPSIZE = parse(Float32, setting)
            end
            if parameter == "STEPS"
                STEPS = parse(Int32, setting)
            end
            if parameter == "LIST_RADIUS"
                LIST_RADIUS = parse(Float32, setting)
            end
            if parameter == "INTERVAL"
                INTERVAL = parse(Int32, setting)
            end
            if parameter == "OUTPUT_FILE"
                OUTPUT_FILE = setting
            end
        end
    end

    settings = SETTINGS(N, DIM, L, INPUT_FILE, TYPE, THERMO, T, CUTOFF, STEPSIZE, STEPS, LIST_RADIUS, INTERVAL, OUTPUT_FILE)

    return settings
end