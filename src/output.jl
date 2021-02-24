basetype(::Outbreak{T}) where {T <: StateType} = T

"""
   collect_results!(results::Dict{Tuple, Any},
      p::Population, t0::Int32, interval::Integer, weight::Integer)

Collect the results of a single sir run on a dictionary of results based on the
final state of population `p`. It records the number of new infections at
regular intervals starting from `t0`. The total number of records is added to
the results with a given `weight`.
"""
function collect_results!(results::Dict{Tuple, Any},
   o::Outbreak, t0::Int32, interval::Int32, weight::Float64)
   # round to nearest interval from the bottom, eg. if t0 is 16:17:35 and we
   # have an interval equal to 10 minutes then we round t0 to 16:10:00
    t0 = t0 - (t0 % interval)
    for (id, states) in o
        # round to the nearest interval from the top, eg. if infection_time is
        # 16:17:35 and we have an interval equal to 10 minutes then we round
        # the infection_time to 16:20:00
        infection_time = time(infection(o, id))
        infection_time = (infection_time + interval - (infection_time % interval)) - t0

        infection_location = location(infection(o, id))
        infection_source = source(infection(o, id))

        if !haskey(results, (:infection, infection_time, infection_source, infection_location))
            results[(:infection, infection_time, infection_source, infection_location)] = 1 * weight
        else
            results[(:infection, infection_time, infection_source, infection_location)] += 1 * weight
        end

        recovery_time = time(state(o, id))
        recovery_time = (recovery_time + interval - (recovery_time % interval)) - t0

        if !haskey(results, (:recovery, recovery_time))
           results[(:recovery, recovery_time)] = 1 * weight
        else
            results[(:recovery, recovery_time)] += 1 * weight
        end
    end
    return results
end

"""
    save(path::AbstractString, results::Dict{Tuple, Any}, params::Dict{String, Any}, N::Int)

Save the simulation results as an HDF5 file.
"""
function save(filename, results::Dict{Tuple,Any}, params::Dict{String,Any}, N::Int)

    infection_time = Vector{Int32}()
    infection_source = Vector{Int32}()
    infection_location = Vector{Int32}()
    infection_value = Vector{Float64}()
    recovery_time = Vector{Int32}()
    recovery_value = Vector{Float64}()
    population_size = -1


    for (ix, (k, v)) in enumerate(results)
        if k[1] == :infection
            push!(infection_time, k[2])
            push!(infection_source, k[2])
            push!(infection_location, k[3])
            push!(infection_value, v)
        elseif k[1] == :recovery
            push!(recovery_time, k[2])
            push!(recovery_value, v)
        elseif k[1] == :population
            population_size = v
        end
    end

    h5open(filename, "w") do file

        write(file, "infection/time", infection_time)
        write(file, "infection/source", infection_source)
        write(file, "infection/location", infection_location)
        write(file, "infection/value", infection_value)

        write(file, "recovery/time", recovery_time)
        write(file, "recovery/value", recovery_value)

        for (k, v) in params
            write(file, "param/$(k)", v)
        end
        write(file, "param/N", N)
        write(file, "param/population_size", population_size)

    end

end

"""
    save(conn::LibPQ.Connection, results::Dict{Tuple, Any}, params::Dict{String, Any}, N::Int)

Save the simulation results as an HDF5 file.
"""
function save(conn::LibPQ.Connection, results::Dict{Tuple,Any}, params::Dict{String,Any}, N::Int)

    params_query = LibPQ.prepare(conn, "SELECT key FROM simulation.params WHERE params = \$1;")
    params_result = LibPQ.execute(params_query, (JSON.json(params),))
    params_key = -1

    if length(params_result) > 0
        params_key = getindex(params_result, 1, 1)
        @info "Deleting previous results, if any."
        LibPQ.execute(conn, "DELETE FROM simulation.results WHERE params_key = \$1;", (params_key,))
    else
        LibPQ.load!(
            (params = [JSON.json(params)], N=[N], population_size=[population_size]),
            conn,
            """
            INSERT INTO simulation.params 
            (params, N, population_size, status) 
            VALUES (\$1, \$2, \$3, 'ongoing');
            """,
        )
        params_result = LibPQ.execute(params_query, (JSON.json(params),))
        params_key = getindex(params_result, 1, 1)
    end

    res_metric = Array{String,1}()
    res_value = Array{Float64,1}()
    res_location = Array{Union{Int32,Missing},1}()
    res_source = Array{Union{Int32,Missing},1}()
    res_elapsed = Array{Int32,1}()
    population_size = -1

    for (ix, (k, v)) in enumerate(results)
        if k[1] == :infection
            push!(res_metric, "infection")
            push!(res_value, v)
            if k[3] == typemax(Int32)
                push!(res_source, missing)
            else
                push!(res_source, k[3])
            end
            push!(res_location, k[4])
            push!(res_elapsed, k[2])
        elseif k[1] == :recovery
            push!(res_metric, "recovery")
            push!(res_value, v)
            push!(res_source, missing)
            push!(res_location, missing)
            push!(res_elapsed, k[2])
        elseif k[1] == :population
            population_size = v
        end
    end

    LibPQ.execute(conn, "BEGIN;")
    try
        LibPQ.load!(
            (
                metric = res_metric,
                value = res_value,
                location = res_location,
                source = res_source,
                elapsed = res_elapsed
            ),
            conn,
            """
            INSERT INTO simulation.results 
            (params_key, metric, value, location, source, elapsed) 
            VALUES ($(params_key), \$1, \$2, \$3, \$4, \$5 * INTERVAL '1 second');
            """,
        )
    finally
        LibPQ.execute(conn, "COMMIT;")
    end

end
