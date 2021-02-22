module Tseir

    using Base
    using HDF5
    using Random
    using StatsBase
    using ProgressBars
    using DataStructures:AbstractMutableHeap
    using Distributions:Sampleable

    import LibPQ
    import Core.Intrinsics.bitcast

    export Individual,
        Population,
        PopulationHeap,
        Model,
        State,
        StateType,

        # methods for Individual
        reset!,
        advance!,
        infect!,
        state,
        infection,
        can_infect,
        contacts,
        transitions,

        # methods for Population
        lazy_init_population,
        eager_init_population,

        # methods for model
        set_rng!,
        set_interevent_distribution!,
        add_path!,
        interevent_distribution,
        advance,

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

        # models
        # SI
        SIR
        # SEIR


    include("./interval.jl")
    include("./state.jl")
    include("./model.jl")
    include("./individual.jl")
    include("./population.jl")
    include("./simulator.jl")
    include("./output.jl")
    include("./sir.jl")

end
