basetype(::Type{<:Individual{T}}) where {T <: StateType} = T

"""
   collect_results!(results::Dict{Tuple{Int32, Int32,Int32},Float64},
      p::Population, t0::Int32, interval::Integer, weight::Integer)

Collect the results of a single sir run on a dictionary of results based on the
final state of population `p`. It records the number of new infections at
regular intervals starting from `t0`. The total number of records is added to
the results with a given `weight`.
"""
function collect_results!(results::Dict{Symbol,Dict{K,Float64} where K},
   p::Population, t0::Int32, interval::Int32, weight::Float64)
   # round to nearest interval from the bottom, eg. if t0 is 16:17:35 and we
   # have an interval equal to 10 minutes then we round t0 to 16:10:00
    t0 = t0 - (t0 % interval)
    for i in p
        if state(i) == typemax(basetype(i))
         # round to the nearest interval from the top, eg. if infection_time is
         # 16:17:35 and we have an interval equal to 10 minutes then we round
         # the infection_time to 16:20:00
            infection_time = time(infection(i))
            infection_time = (infection_time + interval - (infection_time % interval)) - t0

            infection_location = location(infection(i))
            infection_source = source(infection(i))

            if !haskey(results[:infection], (infection_time, infection_source, infection_location))
                results[:infection][(infection_time, infection_source, infection_location)] = 1 * weight
            else
                results[:infection][(infection_time, infection_source, infection_location)] += 1 * weight
            end

            recovery_time = time(state(i))
            recovery_time = (recovery_time + interval - (recovery_time % interval)) - t0

            if !haskey(results[:recovery], recovery_time)
                results[:recovery][recovery_time] = 1 * weight
            else
                results[:recovery][recovery_time] += 1 * weight
            end

        end
    end
    return results
end

"""
   sir_save(path::AbstractString, results::Dict{Tuple{Int32, Int32, Int32}, Float64},
      transmission_probability, recovery_rate, N::Integer, seed::Integer)

Save the simulation results as an HDF5 file.
"""
function save(filename, results::Dict{Symbol,Dict{K,Float64} where K},
   p::Population, m::Model, N::Integer, seed::Integer)

    infection_time = Array{Int32,1}(undef, length(results[:infection]))
    infection_source = Array{Int32,1}(undef, length(results[:infection]))
    infection_location = Array{Int32,1}(undef, length(results[:infection]))
    infection_amount = Array{Float64,1}(undef, length(results[:infection]))

    for (ix, (k, v)) in enumerate(results[:infection])
        infection_time[ix], infection_source[ix], infection_location[ix] = k
        infection_amount[ix] = v
    end

    recovery_time = Array{Int32,1}(undef, length(results[:recovery]))
    recovery_amount = Array{Float64,1}(undef, length(results[:recovery]))

    for (ix, (k, v)) in enumerate(results[:recovery])
        recovery_time[ix] = k
        recovery_amount[ix] = v
    end

    h5open(filename, "w") do file

        write(file, "infection/time", infection_time)
        write(file, "infection/source", infection_source)
        write(file, "infection/location", infection_location)
        write(file, "infection/amount", infection_amount)

        write(file, "recovery/time", recovery_time)
        write(file, "recovery/amount", recovery_amount)

        write(file, "param/transmission_probability", transmission_probability)
        write(file, "param/recovery_rate", recovery_rate)
        write(file, "param/iter", N)
        write(file, "param/seed", seed)
        write(file, "param/population_size", length(p))

    end

end

