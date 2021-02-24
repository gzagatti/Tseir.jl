using Test
using Tseir
using Random
using DataStructures:extract_all!
using Distributions:DiscreteUniform

r = Random.GLOBAL_RNG

function simple_pop()
    p = Population([Individual(id) for id in 1:2])
    p[1].transition_list = [
        Interval(50, 60, typemax(Int32), 1),
        Interval(70, 80, typemax(Int32), 2),
        Interval(85, 99, typemax(Int32), 1),
    ]
    p[2].transition_list = [
        Interval(40, 55, typemax(Int32), 1),
        Interval(75, 95, typemax(Int32), 2),
    ]
    contact_list = [
        (1, 2, 50, 55),
        (1, 2, 75, 77),
        (1, 2, 78, 80),
        (2, 1, 50, 55),
        (2, 1, 75, 77),
        (2, 1, 78, 80),
    ]
    for (i, j, interval_start, interval_end) in contact_list
        Tseir.push_contact!(p[i], p[j].id, interval_start, interval_end)
    end
    return p
end

@testset "Simulator tests" begin

    @testset "transmission_source" begin

        i = Individual(1)
        i.transition_list = [
            Interval(10, 15, typemax(Int32), 1),
            Interval(17, 20, typemax(Int32), 2),
            Interval(32, 41, typemax(Int32), 9),
            Interval(45, 50, typemax(Int32), 2)
        ]
        o = Outbreak{StateSIR}()
        infect!(o, i, Int32(15), Int32(10))
        @test location(infection(o, i)) == 1
        @test migration(infection(o, i)) == 15

        s = Tseir.transmission_source(Int32(11), o, i)
        @test s == 10

        s = Tseir.transmission_source(Int32(18), o, i)
        @test s == 1

        s = Tseir.transmission_source(Int32(47), o, i)
        @test s == 9

    end

    @testset "infect_others!" begin

        p = simple_pop()
        o = Outbreak{StateSIR}()

        h = PopulationHeap(Base.By(i -> time(infection(o, i))))
        infect!(o, p[1], Int32(50), Int32(9))
        push!(h, p[1])
        set_rng!(SIR, r)
        set_interevent_distribution!(SIR, StateSIR(:S), DiscreteUniform(1, 1))
        set_interevent_distribution!(SIR, StateSIR(:I), DiscreteUniform(2, 2))
        Tseir.infect_others!(h, p, SIR, o)
        @test length(h) == 1
        @test type(state(o, p[1])) == StateSIR(:R)
        @test time(state(o, p[1])) == 52
        @test type(state(o, p[2])) == StateSIR(:I)
        @test time(infection(o, p[2])) == 51
        @test source(infection(o, p[2])) == 9
        assign_event_location!(infection(o, p[2]), transitions(p[2]))
        @test location(infection(o, p[2])) == 1

        extract_all!(h)
        empty!(o)
        infect!(o, p[1], Int32(40), Int32(9))
        push!(h, p[1])
        set_rng!(SIR, r)
        set_interevent_distribution!(SIR, StateSIR(:S), DiscreteUniform(1, 1))
        set_interevent_distribution!(SIR, StateSIR(:I), DiscreteUniform(2, 2))
        Tseir.infect_others!(h, p, SIR, o)
        @test type(state(o, p[1])) == StateSIR(:R)
        @test time(state(o, p[1])) == 42
        @test length(h) == 0
        @test type(state(o, p[2])) == StateSIR(:S)
        @test time(infection(o, p[2])) == typemax(Int32)
        @test location(infection(o, p[2])) == typemax(Int32)
        @test source(infection(o, p[2])) == typemax(Int32)

        extract_all!(h)
        empty!(o)
        infect!(o, p[1], Int32(50), Int32(9))
        push!(h, p[1])
        set_rng!(SIR, r)
        set_interevent_distribution!(SIR, StateSIR(:S), DiscreteUniform(1, 1))
        set_interevent_distribution!(SIR, StateSIR(:I), DiscreteUniform(0, 0))
        Tseir.infect_others!(h, p, SIR, o)
        @test type(state(o, p[1])) == StateSIR(:R)
        @test time(state(o, p[1])) == 50
        @test length(h) == 0
        @test type(state(o, p[2])) == StateSIR(:S)
        @test time(infection(o, p[2])) == typemax(Int32)
        @test location(infection(o, p[2])) == typemax(Int32)
        @test source(infection(o, p[2])) == typemax(Int32)

        extract_all!(h)
        empty!(o)
        infect!(o, p[1], Int32(50), Int32(9))
        infect!(o, p[2], Int32(60), Int32(2))
        push!(h, p[1])
        set_rng!(SIR, r)
        set_interevent_distribution!(SIR, StateSIR(:S), DiscreteUniform(1, 1))
        set_interevent_distribution!(SIR, StateSIR(:I), DiscreteUniform(2, 2))
        Tseir.infect_others!(h, p, SIR, o)
        @test type(state(o, p[1])) == StateSIR(:R)
        @test time(state(o, p[1])) == 52
        @test length(h) == 1
        @test first(h) == p[2]
        @test type(state(o, p[2])) == StateSIR(:I)
        @test time(state(o, p[2])) == 51
        @test source(state(o, p[2])) == 9
        assign_event_location!(infection(o, p[2]), transitions(p[2]))
        @test location(infection(o, p[2])) == 1

    end

    @testset "run" begin
        p = simple_pop()
        set_rng!(SIR, r)
        set_interevent_distribution!(SIR, StateSIR(:S), DiscreteUniform(1, 1))
        set_interevent_distribution!(SIR, StateSIR(:I), DiscreteUniform(2, 2))
        o = Tseir.run(p, SIR, p[1], Int32(50))
        @test type(state(o, p[1])) == StateSIR(:R)
        @test type(state(o, p[2])) == StateSIR(:R)
    end

    @testset "sweep" begin
        p = simple_pop()
        set_rng!(SIR, MersenneTwister(100))
        set_interevent_distribution!(SIR, StateSIR(:S), DiscreteUniform(1, 1))
        set_interevent_distribution!(SIR, StateSIR(:I), DiscreteUniform(2, 2))
        results = Tseir.sweep(p, SIR, 50, 55, 2, 1)
        @test results[:infection] == Dict((1, typemax(Int32), 1) => 1.0, (2, typemax(Int32), 1) => 1.0)
        @test results[:recovery] == Dict(4 => 1.0, 3 => 1.0)
    end

end
