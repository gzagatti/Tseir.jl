if !("$(@__DIR__)/../src"    in LOAD_PATH)
    pushfirst!(LOAD_PATH, "$(@__DIR__)/../src")
end

using Tseir
using JSON
using Dates
using LibPQ
using DotEnv
using Printf
using Memento
using ArgParse
using Distributions:Exponential

using RandomNumbers.PCG

basetype(::Model{T}) where {T <: StateType} = T

function main(args)

    params = Dict{String,Any}()

    if args["model"] == "si"
        m = SI
    elseif args["model"] == "sir"
        m = SIR
    elseif args["model"] == "seir"
        m = SEIR
    end

    β = args["beta"]
    set_interevent_distribution!(m, :S, Exponential(1 / β))
    params["beta"] = @sprintf "%0.4e" β

    θ = args["theta"]
    if args["model"] in ["seir"]
        set_interevent_distribution!(m, :E, Exponential(1 / θ))
        params["theta"] = @sprintf "%0.4e" θ
    end

    γ = args["gamma"]
    if args["model"] in ["sir", "seir"]
        set_interevent_distribution!(m, :I, Exponential(1 / γ))
        params["gamma"] = @sprintf "%0.4e" γ
    end

    i0 = args["i0"]
    params["i0"] = "$(i0)"

    epidemics_start = args["start"]
    params["start"] = Dates.format(unix2datetime(epidemics_start), "yyyy-mm-ddTH:M:S")

    epidemics_end = args["end"]
    params["end"] = Dates.format(unix2datetime(epidemics_end), "yyyy-mm-ddTH:M:S")

    T = args["T"]
    params["T"] = T

    N = args["N"]
    params["N"] = T

    save_interval = args["save-interval"]
    params["save_interval"] = save_interval

    seed = args["seed"]
    set_rng!(m, PCG.PCGStateOneseq(seed))
    params["seed"] = seed

    output_root = (
        args["output"] *
        "/$(Dates.format(unix2datetime(args["start"]), "yyyy-mm-dd"))_"    *
        "$(Dates.format(unix2datetime(args["end"]), "yyyy-mm-dd"))"
    )
    mkpath(output_root)

    output_path = args["model"]
    if !isnothing(i0)
        output_path = output_path * "_$(i0)"
    end
    output_path = output_root * "/" * output_path * ".h5"

    eager = args["eager"]

    @info """Simulation parameters:
        env_file: $(args["env-file"])
        db: $(ENV["WIFI_CONN"])
        model: $(m)
        β: $(@sprintf "%0.4e" β)
        θ: $(@sprintf "%0.4e" θ)
        γ: $(@sprintf "%0.4e" γ)
        i0: $(i0)
        epidemics_start: $(Dates.format(unix2datetime(epidemics_start), "yyyy-mm-ddTH:M:S"))
        epidemics_end: $(Dates.format(unix2datetime(epidemics_end), "yyyy-mm-ddTHH:MM:SS"))
        T: $(Dates.format(unix2datetime(T), "yyyy-mm-ddTH:M:S"))
        N: $(@sprintf "%d" N)
        save_interval: $(@sprintf "%d" save_interval)
        seed: $(seed)
        output_path: $(output_path)
        eager: $(args["eager"])
    """

    conn = LibPQ.Connection(ENV["WIFI_CONN"])

    params_query = LibPQ.prepare(conn, "SELECT key FROM simulation.params WHERE params = \$1;")
    params_result = LibPQ.execute(params_query, (JSON.json(params),))
    if length(params_result) > 0
        @info "Simulation with given parameters has already been executed; nothing to run."
        return
    end

    transition_stmt = prepare(
        conn, """
        SELECT
            building_key,
            GREATEST($(epidemics_start), EXTRACT(EPOCH from arrival_time)::INT) AS interval_start,
            LEAST($(T), EXTRACT(EPOCH FROM departure_time)::INT) AS interval_end
        FROM views.bdg_transition
        WHERE userid_key = \$1
        AND departure_time >= '$(Dates.format(unix2datetime(args["start"]), "yyyy-mm-ddTHH:MM:SS"))'
        AND arrival_time < '$(Dates.format(unix2datetime(T), "yyyy-mm-ddTHH:MM:SS"))'
        ORDER BY arrival_time
    """)

    if args["eager"]
        contact_stmt = prepare(
            conn, """
            SELECT
                userid_key,
                userid_key_other,
                GREATEST($(epidemics_start), EXTRACT(EPOCH from overlap_start)::INT),
                LEAST($(T), EXTRACT(EPOCH FROM overlap_end)::INT)
            FROM views.contact_list
            WHERE overlap_end >= '$(Dates.format(unix2datetime(epidemics_start), "yyyy-mm-ddTHH:MM:SS"))'
            AND overlap_start < '$(Dates.format(unix2datetime(T), "yyyy-mm-ddTHH:MM:SS"))'
        """)
        @time p = Tseir.eager_init_population(basetype(m), transition_stmt, contact_stmt)
    else
        id_stmt = prepare(
            conn, """
            SELECT DISTINCT userid_key
            FROM views.bdg_transition
            WHERE departure_time >= '$(Dates.format(unix2datetime(epidemics_start), "yyyy-mm-ddTHH:MM:SS"))'
            AND arrival_time < '$(Dates.format(unix2datetime(T), "yyyy-mm-ddTHH:MM:SS"))'
        """)
        contact_stmt = prepare(
            conn, """
            SELECT
                userid_key_other,
                GREATEST($(epidemics_start), EXTRACT(EPOCH from overlap_start)::INT) AS interval_start,
                LEAST($(T), EXTRACT(EPOCH FROM overlap_end)::INT) AS interval_end
            FROM views.contact_list
            WHERE overlap_end >= '$(Dates.format(unix2datetime(epidemics_start), "yyyy-mm-ddTHH:MM:SS"))'
            AND overlap_start < '$(Dates.format(unix2datetime(T), "yyyy-mm-ddTHH:MM:SS"))'
            AND userid_key = \$1
            ORDER BY overlap_start
        """)
        @time p = Tseir.lazy_init_population(basetype(m), id_stmt, transition_stmt, contact_stmt)
    end


    @info "Population size: $(length(p))"

    results = Dict(
      :infection => Dict{Tuple{Int32,Int32,Int32},Float64}(),
      :recovery => Dict{Int32,Float64}()
    )

    if !isnothing(i0)
        if i0 in p
            i0 = p[i0]
        else
            @info "Patient $(i0) not present in the population given parameters; nothing to run."
            save(output_path, results, params)
            @info "Results saved to: $(output_path)"
            @info "Saving final results to database."
            save(conn, results, params)
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
        sweep_results = sweep!(p, m, epidemics_start, epidemics_end, n_sweep, save_interval, i0)
        for (k, v) in sweep_results[:infection] sweep_results[:infection][k] = v / n_simulation end
        for (k, v) in sweep_results[:recovery] sweep_results[:recovery][k] = v / n_simulation end
        merge!(+, results[:infection], sweep_results[:infection])
        merge!(+, results[:recovery], sweep_results[:recovery])
        save(output_path, results, params)
        @info "Results saved to: $(output_path)"
    end

    @info "Saving final results to database."
    results[:population] = length(p)
    save(conn, results, params)

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
            help = "epidemics start, formatted as y-m-d H:M:S"
            default = "2020-01-13 00:00:00"
        "--end"
            help = "epidemics end, formatted as y-m-d H:M:S"
            default = "2020-01-24 23:59:59"
        "--T"
            help = "simulation end, formatted as y-m-d H:M:S"
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
            help = "output root"
            default = "$(get(ENV, "HOST_DATA", "/tmp"))/simulations-v5"
        "--eager"
            help = "whether to perform eager or lazy population initialization"
            action = :store_true
    end

    args = parse_args(s)

    if !isnothing(args["env-file"])
        DotEnv.config(args["env-file"])
        ENV["WIFI_CONN"] =
            "postgres://agens:" *
            "$(ENV["AGENS_PW"])@0.0.0.0:"    *
            "$(ENV["AGENS_PORT"])/wifidb"
    end

    args["model"] = lowercase(args["model"])
    args["output"] = "$(get(ENV, "HOST_DATA", "/tmp"))/simulations-v5"
    args["start"] = Int32(datetime2unix(DateTime(args["start"], "y-m-d H:M:S")))
    args["end"] = Int32(datetime2unix(DateTime(args["end"], "y-m-d H:M:S")))
    args["T"] = Int32(datetime2unix(DateTime(args["T"], "y-m-d H:M:S")))

    main(args)

end

cli()
