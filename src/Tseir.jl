module Tseir

    include("./agents.jl")
    include("./simulator.jl")

    export Individual, Population, PopulationHeap
    export reset!, update!, transitions, contacts
    export lazy_init_population, eager_init_population
    export sir!, sir_sweep!, sir_save

end
