module Tseir

    using Base
    using HDF5
    using Enums
    using Random
    using StatsBase
    using ProgressBars
    using DataStructures: AbstractMutableHeap

    import LibPQ
    import Core.Intrinsics.bitcast

    include("./individuals.jl")
    include("./models.jl")
    include("./output.jl")
    include("./population.jl")
    include("./simulator.jl")
    include("./sir.jl")
    include("./states.jl")

    export 

        # structures
        Individual,
        Population,
        PopulationHeap,
        Model,
        State,
        StateType

        # models
        # SI
        SIR
        # SEIR

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

        # methods for StateType, State
        @states,

        type,
        time,
        location,
        source,
        migration

end
