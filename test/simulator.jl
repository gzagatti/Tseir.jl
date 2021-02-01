if !("$(@__DIR__)/../src" in LOAD_PATH)
    pushfirst!(LOAD_PATH, "$(@__DIR__)/../src")
end

using Test
using Tseir
using DataStructures: extract_all!
using Distributions: DiscreteUniform

r = Random.GLOBAL_RNG

@testset "Simulator tests" begin

    @testset "transmission" begin

        contact_events = [
            Tseir.ContactEvent(1, 10, 15, 5),
            Tseir.ContactEvent(2, 17, 20, 8),
            Tseir.ContactEvent(9, 32, 41, 17),
            Tseir.ContactEvent(2, 45, 50, 22)
        ]

        t = Tseir.transmission(contact_events, Int32(5), DiscreteUniform(1, 1), r) 
        @test t == 11

        t = Tseir.transmission(contact_events, Int32(18), DiscreteUniform(2, 2), r)
        @test t == 20

        t = Tseir.transmission(contact_events, Int32(19), DiscreteUniform(11, 11), r)
        @test t == 46

        t = Tseir.transmission(contact_events, Int32(15), DiscreteUniform(2, 2), r)
        @test t == 19

        t = Tseir.transmission(contact_events, Int32(52), DiscreteUniform(1, 1), r)
        @test t == nothing

        t = Tseir.transmission(contact_events, Int32(30), DiscreteUniform(50, 50), r)
        @test t == nothing

    end

    @testset "transmission_source" begin

    i = Individual(1)
    i.migration_time = Int32(15)
    i.infection_source = 10
    i.transition_list = [
        Tseir.TransitionEvent(1, 10, 15),
        Tseir.TransitionEvent(2, 17, 20),
        Tseir.TransitionEvent(9, 32, 41),
        Tseir.TransitionEvent(2, 45, 50)
    ]

    s = Tseir.transmission_source(Int32(11), i)
    @test s == 10

    s = Tseir.transmission_source(Int32(18), i)
    @test s == 1

    s = Tseir.transmission_source(Int32(47), i)
    @test s == 9

    end

    @testset "infect_others!" begin
        p = Population([Individual(id) for id in 1:2])
        p[1].transition_list = [
           Tseir.TransitionEvent(1, 50, 60),
           Tseir.TransitionEvent(2, 70, 80),
           Tseir.TransitionEvent(1, 85, 99),
        ]
        p[2].transition_list = [
           Tseir.TransitionEvent(1, 40, 55),
           Tseir.TransitionEvent(2, 75, 95),
        ]
        contact_list = [
            (1, 2, 1, 50, 55),
            (1, 2, 4, 75, 77),
            (1, 2, 5, 78, 80),
            (2, 1, 1, 50, 55),
            (2, 1, 4, 75, 77),
            (2, 1, 5, 78, 80),
        ]
        for (i, j, coord, event_start, event_end) in contact_list
            Tseir.push_contact!(p[i], p[j].id, coord, event_start, event_end)
        end

        h = PopulationHeap(Base.By(i -> i.infection_time))
        p[1].infection_time = 50
        p[1].infection_source = 9
        push!(h, p[1])
        Tseir.infect_others!(h, p, DiscreteUniform(1, 1), DiscreteUniform(2, 2), r)
        @test p[1].status == Tseir.R
        @test p[1].recovery_time == 52 
        @test length(h) == 1
        @test first(h) == p[2]
        @test p[2].status == Tseir.I
        @test p[2].infection_time == 51
        @test p[2].infection_location == 1
        @test p[2].infection_source == 9

        extract_all!(h)
        reset!(p)
        p[1].infection_time = 40
        p[1].infection_source = 9
        push!(h, p[1])
        Tseir.infect_others!(h, p, DiscreteUniform(1, 1), DiscreteUniform(2, 2), r)
        @test p[1].status == Tseir.R
        @test p[1].recovery_time == 42 
        @test length(h) == 0
        @test p[2].status == Tseir.S
        @test p[2].infection_time == typemax(Int32)
        @test p[2].infection_location == typemax(Int32)
        @test p[2].infection_source == typemax(Int32)

        extract_all!(h)
        reset!(p)
        p[1].infection_time = 50
        p[1].infection_source = 9
        push!(h, p[1])
        Tseir.infect_others!(h, p, DiscreteUniform(1, 1), DiscreteUniform(0, 0), r)
        @test p[1].status == Tseir.R
        @test p[1].recovery_time == 50
        @test length(h) == 0
        @test p[2].status == Tseir.S
        @test p[2].infection_time == typemax(Int32)
        @test p[2].infection_location == typemax(Int32)
        @test p[2].infection_source == typemax(Int32)

        h = PopulationHeap(Base.By(i -> i.infection_time))
        p[1].infection_time = 50
        p[1].infection_source = 9
        p[2].infection_time = 60
        push!(h, p[1])
        Tseir.infect_others!(h, p, DiscreteUniform(1, 1), DiscreteUniform(2, 2), r)
        @test p[1].status == Tseir.R
        @test p[1].recovery_time == 52 
        @test length(h) == 1
        @test first(h) == p[2]
        @test p[2].status == Tseir.I
        @test p[2].infection_time == 51
        @test p[2].infection_location == 1
        @test p[2].infection_source == 9

    end

end
