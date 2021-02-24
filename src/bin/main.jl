using Tseir
using Dates
using DotEnv
using Printf
using Memento
using ArgParse
using Distributions:Exponential

import LibPQ
import JSON
import RandomNumbers.PCG.PCGStateOneseq

function simulate(params::Dict{String, Any}, N::Int, conn::LibPQ.Connection, eager::Bool,
    p::Union{Population, Nothing}=nothing, output_path::Union{String, Nothing}=nothing)

    results = Dict{Tuple, Any}()
    params_key = -1

    if params["model"] == "si"
        m = SI
    elseif params["model"] == "sir"
        m = SIR
    elseif params["model"] == "seir"
        m = SEIR
    end

    set_interevent_distribution!(m, :S, Exponential(1 / parse(Float64, params["beta"])))

    if params["model"] in ["seir"]
        set_interevent_distribution!(m, :E, Exponential(1 / parse(Float64, params["theta"])))
    end

    if params["model"] in ["sir", "seir"]
        set_interevent_distribution!(m, :I, Exponential(1 / parse(Float64, params["gamma"])))
    end

    i0 = params["i0"] == "nothing" ? nothing : parse(Int32, params["i0"])

    epidemics_start = Int32(datetime2unix(DateTime(params["start"], "yyyy-mm-dd HH:MM:SS")))
    epidemics_end = Int32(datetime2unix(DateTime(params["end"], "yyyy-mm-dd HH:MM:SS")))
    T = Int32(datetime2unix(DateTime(params["T"], "yyyy-mm-dd HH:MM:SS")))

    save_interval = params["save_interval"]

    set_rng!(m, PCGStateOneseq(params["seed"]))

    msg = """Simulation parameters:
        conn: $(conn)
        model: $(m)
        β: $(params["beta"])
    """

    if haskey(params, "theta")
        msg *= "θ: $(params["theta"])\n"
    end

    if haskey(params, "gamma")
        msg *= "γ: $(params["gamma"])\n"
    end

    msg *= """i0: $(i0)
        epidemics_start: $(params["start"])
        epidemics_end: $(params["end"])
        T: $(params["T"])
        N: $(N)
        save_interval: $(params["save_interval"])
        seed: $(params["seed"])
        output_path: $(output_path)
        eager: $(eager)
    """

    params_query = LibPQ.prepare(conn, "SELECT key, n, status FROM simulation.params WHERE params = \$1;")
    params_result = LibPQ.execute(params_query, (JSON.json(params),))

    if length(params_result) > 0
        params_key = getindex(params_result, 1, 1)
        if (N <= getindex(params_result, 1, 2)) && (getindex(params_result, 1, 3) != "error")
            @info "Simulation with given parameters has already been executed; nothing to run."
            return
        end
    end

    try
        if length(params_result) > 0
            LibPQ.execute(
                conn,
                """
                UPDATE simulation.params SET
                status = 'ongoing'
                WHERE key = \$1
                """,
                (params_key,),
            )
        else
            LibPQ.load!(
                (params = [JSON.json(params)], N=[N], population_size=[missing]),
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

        if isnothing(p)
            p = cli_init_population(conn, epidemics_start, T, eager)
        end

        @info "Total population size: $(length(p))"
        results[(:population,)] = length(p)

        if !isnothing(i0)
            if i0 in p
                i0 = p[i0]
            else
                @info "Patient $(i0) not present in the population given parameters; nothing to run."
                if isnothing(output_path)
                    @info "Saving results to database."
                    save(conn, results, params, N)
                else
                    @info "Saving results to: $(output_path)"
                    save(output_path, results, params, N)
                    @info "Updating status to done."
                    LibPQ.execute(
                        conn,
                        """
                        UPDATE simulation.params SET
                        n = \$1, population_size = \$2, status = 'done'
                        WHERE key = \$3
                        """,
                        (N, length(p), params_key),
                    )
                end
                return
            end
        end

        n_sweep = 100
        n_simulation = ceil(Int, N / n_sweep)
        start = now()
        for ix in 1:n_simulation
            if (ix == n_simulation) && (N % n_sweep != 0)
                n_sweep = N % n_sweep
            end
            @info "Simulation $ix / $(n_simulation) ($(canonicalize(Dates.CompoundPeriod(now() - start))) elapsed)"
            sweep_results = sweep(p, m, epidemics_start, epidemics_end, n_sweep, save_interval, i0)
            for (k, v) in sweep_results sweep_results[k] = v / (N/n_sweep) end
            merge!(+, results, sweep_results)
            if isnothing(output_path)
                @info "Saving final results to database."
                save(conn, results, params, N)
            else
                @info "Saving results to: $(output_path)"
                save(output_path, results, params, N)
            end
        end

        @info "Updating status to done."
        LibPQ.execute(
            conn,
            """
            UPDATE simulation.params SET
            n = \$1, population_size = \$2, status = 'done'
            WHERE key = \$3
            """,
            (N, length(p), params_key),
        )

    catch
        if params_key > -1
            LibPQ.execute(
                conn,
                """
                UPDATE simulation.params SET
                status = 'error'
                WHERE key = \$1
                """,
                (params_key,),
            )
        end
    end

end

function cli_init_population(conn::LibPQ.Connection, epidemics_start::Int32, T::Int32, eager::Bool)

        transition_stmt = LibPQ.prepare(
            conn, """
            SELECT
                building_key,
                GREATEST($(epidemics_start), EXTRACT(EPOCH from arrival_time)::INT) AS interval_start,
                LEAST($(T), EXTRACT(EPOCH FROM departure_time)::INT) AS interval_end
            FROM views.bdg_transition
            WHERE userid_key = \$1
            AND departure_time >= '$(Dates.format(unix2datetime(epidemics_start), "yyyy-mm-dd HH:MM:SS"))'
            AND arrival_time < '$(Dates.format(unix2datetime(T), "yyyy-mm-dd HH:MM:SS"))'
            ORDER BY arrival_time
        """)

        if eager
            contact_stmt = LibPQ.prepare(
                conn, """
                SELECT
                    userid_key,
                    userid_key_other,
                    GREATEST($(epidemics_start), EXTRACT(EPOCH from overlap_start)::INT),
                    LEAST($(T), EXTRACT(EPOCH FROM overlap_end)::INT)
                FROM views.contact_list
                WHERE overlap_end >= '$(Dates.format(unix2datetime(epidemics_start), "yyyy-mm-dd HH:MM:SS"))'
                AND overlap_start < '$(Dates.format(unix2datetime(T), "yyyy-mm-dd HH:MM:SS"))'
            """)
            @time p = Tseir.eager_init_population(transition_stmt, contact_stmt)
        else
            id_stmt = LibPQ.prepare(
                conn, """
                SELECT DISTINCT userid_key
                FROM views.bdg_transition
                WHERE departure_time >= '$(Dates.format(unix2datetime(epidemics_start), "yyyy-mm-ddTHH:MM:SS"))'
                AND arrival_time < '$(Dates.format(unix2datetime(T), "yyyy-mm-dd HH:MM:SS"))'
            """)
            contact_stmt = LibPQ.prepare(
                conn, """
                SELECT
                    userid_key_other,
                    GREATEST($(epidemics_start), EXTRACT(EPOCH from overlap_start)::INT) AS interval_start,
                    LEAST($(T), EXTRACT(EPOCH FROM overlap_end)::INT) AS interval_end
                FROM views.contact_list
                WHERE overlap_end >= '$(Dates.format(unix2datetime(epidemics_start), "yyyy-mm-dd HH:MM:SS"))'
                AND overlap_start < '$(Dates.format(unix2datetime(T), "yyyy-mm-dd HH:MM:SS"))'
                AND userid_key = \$1
                ORDER BY overlap_start
            """)
            @time p = Tseir.lazy_init_population(id_stmt, transition_stmt, contact_stmt)
        end


        @info "Population size: $(length(p))"

    return p

end

function cli()

    s = ArgParseSettings()
    @add_arg_table s begin
        "model"
            help = "model, one of [si, sir, seir]"
            range_tester = x -> lowercase(x) in ["si", "sir", "seir"]
            required = true
            arg_type = String
        "--env-file", "-e"
            help = "environment file with database connection settings and path to host data."
        "--beta", "-b"
            help = "transmission rate (per time units of the model)"
            arg_type = Float64
            default = 1 / (60 * 60 * 24 * 3)
        "--theta", "-t"
            help = "incubation rate (per time units of the model)"
            arg_type = Float64
            default = 1 / (60 * 60 * 24 * 3)
        "--gamma", "-g"
            help = "recovery rate (per time units of the model)"
            arg_type = Float64
            default = 1 / (60 * 60 * 24 * 14)
        "--i0"
            help = "patient zero"
            default = nothing
            arg_type = Int32
        "--start"
            help = "epidemics start, formatted as yyyy-mm-dd HH:MM:SS"
            default = "2020-01-13 00:00:00"
        "--end"
            help = "epidemics end, formatted as yyyy-mm-dd HH:MM:SS"
            default = "2020-01-24 23:59:59"
        "--T"
            help = "simulation end, formatted as yyyy-mm-dd HH:MM:SS"
            default = "2020-04-06 23:59:59"
        "--N"
            help = "number of simulations"
            default = 1000
            arg_type = Int
        "--save-interval", "-i"
            help = "save interval (default every 10 minutes)"
            default = 15 * 60
            arg_type = Number
        "--seed", "-s"
            help = "random seed"
            default = 5921
            arg_type = Int
        "--output", "-o"
            help = "output path, if specified save results to file instead of database."
        "--eager"
            help = "whether to perform eager or lazy population initialization"
            action = :store_true
    end

    args = parse_args(s)

    args["start"] = Int32(datetime2unix(DateTime(args["start"], "yyyy-mm-dd HH:MM:SS")))
    args["end"] = Int32(datetime2unix(DateTime(args["end"], "yyyy-mm-dd HH:MM:SS")))
    args["T"] = Int32(datetime2unix(DateTime(args["T"], "yyyy-mm-dd HH:MM:SS")))

    params = Dict{String,Any}()

    params["model"] = lowercase(args["model"])
    params["beta"] = @sprintf "%0.8e" args["beta"]
    if params["model"] in ["seir"]
        params["theta"] = @sprintf "%0.8e" args["theta"]
    end
    if args["model"] in ["sir", "seir"]
        params["gamma"] = @sprintf "%0.8e" args["gamma"]
    end
    params["i0"] = "$(args["i0"])"
    params["seed"] = args["seed"]

    params["start"] = Dates.format(unix2datetime(args["start"]), "yyyy-mm-dd HH:MM:SS")
    params["end"] = Dates.format(unix2datetime(args["end"]), "yyyy-mm-dd HH:MM:SS")
    params["T"] = Dates.format(unix2datetime(args["T"]), "yyyy-mm-dd HH:MM:SS")

    params["save_interval"] = args["save-interval"]

    if !isnothing(args["env-file"])
        DotEnv.config(args["env-file"])
        ENV["WIFI_CONN"] =
            "postgres://agens:" *
            "$(ENV["AGENS_PW"])@0.0.0.0:"    *
            "$(ENV["AGENS_PORT"])/wifidb"
    end

    conn = LibPQ.Connection(ENV["WIFI_CONN"])

    if !isnothing(args["output"])
        simulate(params, args["N"], conn, args["eager"], nothing, args["output"])
    end

    simulate(params, args["N"], conn, args["eager"], nothing, nothing)

end

if abspath(PROGRAM_FILE) == @__FILE__
    cli()
end
