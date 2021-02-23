if !(realpath("$(@__DIR__)/../Tseir.jl/src") in LOAD_PATH)
    pushfirst!(LOAD_PATH, realpath("$(@__DIR__)/../Tseir.jl/src"))
end

using Tseir
using Dates
using LibPQ
using DotEnv
using Printf
using Memento
using ArgParse

using RandomNumbers.PCG

function main(args)

    epidemics_start = args["start"]
    epidemics_end = args["end"]
    β = args["beta"]
    γ = args["gamma"]
    N = args["N"]
    seed = args["seed"]

    output_root = (
        args["output"] *
        "/$(Dates.format(unix2datetime(args["start"]), "yyyy-mm-dd"))_"   *
        "$(Dates.format(unix2datetime(args["end"]), "yyyy-mm-dd"))"
    )
    mkpath(output_root)

    output_path = @sprintf "tsir_%0.1e_%0.1e_%0.1e" N β γ
    output_path = output_root * "/" * output_path * ".h5"

    @info """Simulation parameters:
        env_file: $(args["env-file"])
        db: $(ENV["WIFI_CONN"])
        epidemics_start: $(Dates.format(unix2datetime(epidemics_start), "yyyy-mm-ddTH:M:S"))
        epidemics_end: $(Dates.format(unix2datetime(epidemics_end), "yyyy-mm-ddTHH:MM:SS"))
        β: $(@sprintf "%0.4e" β)
        γ: $(@sprintf "%0.4e" γ)
        N: $(@sprintf "%d" N)
        seed: $(seed)
        eager: $(args["eager"])
        output_path: $(output_path)
    """

    conn = LibPQ.Connection(ENV["WIFI_CONN"])

    transition_stmt = prepare(
        conn, """
        SELECT
            bdg_fl_key,
            GREATEST($(args["start"]), EXTRACT(EPOCH from arrival_time)::INT),
            EXTRACT(EPOCH FROM departure_time)::INT
        FROM views.bdg_fl_transition 
        WHERE userid_key = \$1
        AND departure_time >= '$(Dates.format(unix2datetime(args["start"]), "yyyy-mm-ddTHH:MM:SS"))'
    """)

    if args["eager"]
        contact_stmt = prepare(
            conn, """
            SELECT
                userid_key,
                userid_key_other,
                ap_key,
                GREATEST($(args["start"]), EXTRACT(EPOCH from overlap_start)::INT),
                EXTRACT(EPOCH FROM overlap_end)::INT
            FROM views.contact_list
            WHERE overlap_end >= '$(Dates.format(unix2datetime(args["start"]), "yyyy-mm-ddTHH:MM:SS"))'
        """)
        @time p = Tseir.eager_init_population(transition_stmt, contact_stmt)
    else
        id_stmt = prepare(
            conn, """
            SELECT DISTINCT userid_key
            FROM views.bdg_fl_transition 
            WHERE departure_time >= '$(Dates.format(unix2datetime(args["start"]), "yyyy-mm-ddTHH:MM:SS"))'
        """)
        contact_stmt = prepare(
            conn, """
            SELECT
                userid_key_other,
                ap_key,
                GREATEST($(args["start"]), EXTRACT(EPOCH from overlap_start)::INT),
                EXTRACT(EPOCH FROM overlap_end)::INT
            FROM views.contact_list
            WHERE overlap_end >= '$(Dates.format(unix2datetime(args["start"]), "yyyy-mm-ddTHH:MM:SS"))'
            AND userid_key = \$1
        """)
        @time p = Tseir.lazy_init_population(id_stmt, transition_stmt, contact_stmt)
    end


    @info "Population size: $(length(p))"

    r = PCG.PCGStateOneseq(seed)
    results = Dict(
      :infection => Dict{Tuple{Int32,Int32,Int32},Float64}(),
      :recovery => Dict{Int32,Float64}()
    )
    n_sweep = 100
    n_simulation = floor(Int, N / n_sweep)
    start = now()
    for ix in 1:n_simulation
        @info "Simulation $ix / $(n_simulation) ($(canonicalize(Dates.CompoundPeriod(now() - start))) elapsed)"
        sweep_results = sir_sweep!(p, epidemics_start, epidemics_end, β, γ, n_sweep, r)
        for (k, v) in sweep_results[:infection] sweep_results[:infection][k] = v / n_simulation end
        for (k, v) in sweep_results[:recovery] sweep_results[:recovery][k] = v / n_simulation end
        merge!(+, results[:infection], sweep_results[:infection])
        merge!(+, results[:recovery], sweep_results[:recovery])
        sir_save(output_path, results, p, β, γ, ix * ix, seed)
        @info "Results saved to: $(output_path)"
        println()
    end

end

function cli()

    s = ArgParseSettings()
    @add_arg_table s begin
        "--env-file", "-e"
        help = "environment file with database connection settings and path to host data."
        default = realpath("$(@__DIR__)/../../.env")
        "--beta", "-b"
        help = "transmission probability (per time units of the model)"
        arg_type = Float64
        default = 1 / (60 * 60)
        "--gamma", "-g"
        help = "recovery rate (per time units of the model)"
        arg_type = Float64
        default = 1 / (14 * 24 * 60 * 60)
        "--start"
            help = "epidemics start, formatted as y-m-d H:M:S"
            default = "2020-01-13 00:00:00"
        "--end"
            help = "epidemics end, formatted as y-m-d H:M:S"
            default = "2020-01-24 23:59:59"
        "--N"
            help = "number of simulations"
            default = 1000
            arg_type = Int
        "--seed", "-s"
            help = "random seed"
            default = 5921
            arg_type = Int
        "--output", "-o"
            help = "output root"
            default = "$(get(ENV, "HOST_DATA", "/tmp"))/simulations-v4"
        "--eager"
            help = "whether to perform eager or lazy population initialization"
            action = :store_true
    end

    args = parse_args(s)
    DotEnv.config(args["env-file"])

    ENV["WIFI_CONN"] =
        "postgres://agens:" *
        "$(ENV["AGENS_PW"])@0.0.0.0:"   *
        "$(ENV["AGENS_PORT"])/wifidb"

    args["output"] = "$(get(ENV, "HOST_DATA", "/tmp"))/simulations-v4"
    args["start"] = Int32(datetime2unix(DateTime(args["start"], "y-m-d H:M:S")))
    args["end"] = Int32(datetime2unix(DateTime(args["end"], "y-m-d H:M:S")))

    main(args)

end

cli()
