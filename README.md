# Tseir.jl: Temporal SIR

Temporal susceptible -> exposed -> infected -> recovered (Tseir) simulator, based on 
[pholme/tsir](https://github.com/pholme/tsir).

The simulator is developed as a Julia package.

## Similar projects

* *[Pathogen.jl](https://github.com/jangevaare/Pathogen.jl)*: Simulation,
visualization, and inference tools for modelling the spread of infectious
diseases with Julia.
    * keeps a list of transmission and event rates to sample the next event.
    * keeps a transmission matrix between individuals, while we
      compute whether a transmission event takes place when popping from the
      transmission list.
* *[BioSimulator.jl](https://github.com/alanderos91/BioSimulator.jl/tree/master)*: stochastic simulation of complex biology systems
