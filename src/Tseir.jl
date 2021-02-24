module Tseir

    using Base
    using HDF5
    using Random
    using Random: AbstractRNG
    using StatsBase
    using ProgressBars
    using DataStructures:AbstractMutableHeap
    using Distributions:Sampleable

    import JSON
    import LibPQ
    import Core.Intrinsics.bitcast

    export Individual,
        Population,
        PopulationHeap,
        Model,
        Outbreak,
        State,
        StateType,
        Interval,

        # methods for Individual
        contacts,
        transitions,

        # methods for Population
        lazy_init_population,
        eager_init_population,

        # method for PopulationHeap
        update!,

        # methods for Model
        set_rng!,
        set_interevent_distribution!,
        add_path!,
        interevent_distribution,
        advance,

        # methods for Outbreak
        advance!,
        infect!,
        rewind!,
        state,
        previous_state,
        infection,
        susceptible,

        # methods for StateType
        @states,
        exposed,
        infectious,

        # methods for State
        type ,
        time,
        location,
        source,
        migration,
        assign_event_location!,

        # excutables
        run,
        sweep,
        collect_results!,
        save,
        simulate,

        # models
        SI,
        StateSI,
        SIR,
        StateSIR,
        SEIR,
        StateSEIR


    include("./interval.jl")
    include("./state.jl")
    include("./model.jl")
    include("./individual.jl")
    include("./outbreak.jl")
    include("./population.jl")
    include("./simulator.jl")
    include("./output.jl")
    include("./models/si.jl")
    include("./models/sir.jl")
    include("./models/seir.jl")
    include("./bin/main.jl")

end
