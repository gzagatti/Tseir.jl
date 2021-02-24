fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"
quiet = length(ARGS) > 0 && ARGS[1] == "-q"
coverage = length(ARGS) > 0 && ARGS[1] == "-c"
anyerrors = false

using Tseir

if !coverage
    using Test
    my_tests = [
        "individual.jl",
        "population.jl",
        "outbreak.jl",
        "model.jl",
        "simulator.jl"
    ]

    println("Running tests:")
    println()

    for t in my_tests
        try
            include(t)
            println("\033[1m\033[32mPASSED\033[0m: $(t)")
            println()
        catch e
            global anyerrors = true
            println("\033[1m\033[31mFAILED\033[0m: $(t)")
            println()
            if fatalerrors
                rethrow(e)
            elseif !quiet
                showerror(stdout, e, backtrace())
                println()
            end
        end
    end
end

if coverage
    using Coverage
    coverage = process_folder("$(@__DIR__)/../src")
    LCOV.writefile("$(@__DIR__)/coverage-lcov.info", coverage)
    clean_folder("$(@__DIR__)/..")
    coverage = merge_coverage_counts(coverage, LCOV.readfolder(@__DIR__))
    covered_lines, total_lines = get_summary(coverage)
    fraction_covered = round(100 * (covered_lines / total_lines), digits=1)
    println("Test coverage: $(fraction_covered)% ($(covered_lines) / $(total_lines))")
end

