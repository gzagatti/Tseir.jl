if !("$(@__DIR__)/../src"  in LOAD_PATH)
    pushfirst!(LOAD_PATH, "$(@__DIR__)/../src")
end

using Test
using Tseir
using Random
using Distributions:DiscreteUniform

r = Random.GLOBAL_RNG

@testset "Agents tests" begin

    @testset "Individual" begin
        individuals = [Individual(id) for id in 1:2]

        for ix in 1:length(individuals)
    @test individuals[ix].id == ix
    @test individuals[ix].status == Tseir.S
    @test individuals[ix].infection_time == typemax(Int32)
    @test individuals[ix].infection_location == typemax(Int32)
    @test individuals[ix].infection_source == typemax(Int32)
    @test individuals[ix].migration_time == typemax(Int32)
    @test individuals[ix].recovery_time == typemax(Int32)
    @test !isdefined(individuals[ix], :transition_list)
    @test !isdefined(individuals[ix], :contact_list)
end

        transition_list = [
           Tseir.TransitionEvent(1, 50, 60),
           Tseir.TransitionEvent(2, 70, 80),
           Tseir.TransitionEvent(1, 85, 99),
        ]
        Tseir.cat_transition!(individuals[1], transition_list)

        @test length(individuals[1].transition_list) == 3

        transition_list = [
           Tseir.TransitionEvent(1, 40, 55),
           Tseir.TransitionEvent(2, 75, 95),
        ]
        Tseir.cat_transition!(individuals[2], transition_list)

        @test length(individuals[2].transition_list) == 2

        for i in individuals
    for ix in 2:length(i.transition_list)
        @test i.transition_list[ix - 1].event_end < i.transition_list[ix].event_start
        @test i.transition_list[ix - 1].coordset != i.transition_list[ix].coordset
    end
end

        i1 = individuals[1]
        i2 = individuals[2]

        contact_list = [
            (i1, i2.id, 4, 75, 77),
            (i1, i2.id, 1, 50, 55),
            (i1, i2.id, 5, 78, 80),
            (i2, i1.id, 4, 75, 77),
            (i2, i1.id, 5, 78, 80),
            (i2, i1.id, 1, 50, 55),
        ]
        for (i, otherid, coord, event_start, event_end) in contact_list
    Tseir.push_contact!(i, otherid, coord, event_start, event_end)
end

        for i in individuals
    Tseir.sort_contacts!(i)
    for ix in 2:length(i.contact_list)
        @test i.contact_list[ix - 1] <= i.contact_list[ix]
    end
end

        @test haskey(i1.contact_list, i2.id)
        @test haskey(i2.contact_list, i1.id)

        @test i1.contact_list[i2.id][1].coord == 1
        @test i1.contact_list[i2.id][1].event_start == 50
        @test i1.contact_list[i2.id][1].event_end == 55
        @test i1.contact_list[i2.id][1].cum == 5

        @test i1.contact_list[i2.id][2].coord == 4
        @test i1.contact_list[i2.id][2].event_start == 75
        @test i1.contact_list[i2.id][2].event_end == 77
        @test i1.contact_list[i2.id][2].cum == 7

        @test i1.contact_list[i2.id][3].coord == 5
        @test i1.contact_list[i2.id][3].event_start == 78
        @test i1.contact_list[i2.id][3].event_end == 80
        @test i1.contact_list[i2.id][3].cum == 9

        Tseir.infect!(i1, 45, 3)
        @test i1.status == Tseir.I
        @test i1.infection_time == 45
        @test i1.infection_source == 3
        @test i1.infection_location == 1
        @test i1.migration_time == 60

        Tseir.infect!(i1, 50, 4)
        @test i1.status == Tseir.I
        @test i1.infection_time == 50
        @test i1.infection_source == 4
        @test i1.infection_location == 1
        @test i1.migration_time == 60

        Tseir.infect!(i1, 80, 3)
        @test i1.status == Tseir.I
        @test i1.infection_time == 80
        @test i1.infection_source == 3
        @test i1.infection_location == 2
        @test i1.migration_time == 80

        Tseir.infect!(i1, 82, 3)
        @test i1.status == Tseir.I
        @test i1.infection_time == 82
        @test i1.infection_source == 3
        @test i1.infection_location == 1
        @test i1.migration_time == 99

        Tseir.infect!(i1, 105, 6)
        @test i1.status == Tseir.I
        @test i1.infection_time == 105
        @test i1.infection_source == 6
        @test i1.infection_location == typemax(Int32)
        @test i1.migration_time == typemax(Int32)

        reset!(i1)
        @test i1.status == Tseir.S
        @test i1.infection_time == typemax(Int32)
        @test i1.infection_location == typemax(Int32)
        @test i1.infection_source == typemax(Int32)
        @test i1.migration_time == typemax(Int32)
        @test i1.recovery_time == typemax(Int32)

        i1.infection_time = 5
        Tseir.recover!(i1, DiscreteUniform(10, 10), r)
        @test i1.recovery_time == 15
        @test i1.status == Tseir.R

        reset!(i1)

    end

    @testset "Population" begin

        p = Population()

        @test isempty(p) == true
        @test length(p) == 0

        for i in 1:7
    push!(p, Individual(i))
end

        @test isempty(p) == false
        @test length(p) == 7

        @test 1 in p
        @test Individual(1) in p

        i = p[2]
        @test typeof(i) == Individual
        @test i.id == 2

        i = p[Individual(2)]
        @test typeof(i) == Individual
        @test i.id == 2

        push!(p, Individual(1))
        @test length(p) == 7

        i = pop!(p, Individual(1))
        @test i.id == 1
        @test typeof(i) == Individual
        @test length(p) == 6

        i = pop!(p, 2)
        @test i.id == 2
        @test typeof(i) == Individual
        @test length(p) == 5

        q = delete!(p, Individual(3))
        @test q == p
        @test length(p) == 4

        delete!(p, Individual(3))
        @test length(p) == 4

        q = delete!(p, 4)
        @test q == p
        @test length(p) == 3

        delete!(p, 4)
        @test length(p) == 3

        i = pop!(p)
        @test typeof(i) == Individual
        @test length(p) == 2

        for i in p
    i.status = Tseir.I
end

        reset!(p)
        for i in p
    @test i.status == Tseir.S
end

        empty!(p)
        @test isempty(p) == true

    end

    @testset "PopulationHeap" begin

        h = PopulationHeap(Base.By(i -> i.infection_time))
        p = Population()

        @test isempty(h) == true
        @test length(h) == 0

        for id in 1:7
    i = Individual(id)
    i.infection_time = floor(Int32, 100 / id)
    push!(p, i)
    push!(h, i)
end

        @test isempty(h) == false
        @test length(h) == 7

        @test Individual(3) in h
        @test 3 in h

        i = h[Individual(2)]
        @test typeof(i) == Individual
        @test i.id == 2

        i = h[2]
        @test typeof(i) == Individual
        @test i.id == 2
        @test i == p[2]

        i = first(h)
        @test typeof(i) == Individual
        @test i.id == 7
        @test i == p[7]

        p[1].infection_time = 10
        update!(h, p[1])
        @test length(h) == 7
        @test first(h).id == 1
        @test first(h) == p[1]

        p[2].infection_time = 8
        update!(h, 2)
        @test length(h) == 7
        @test first(h).id == 2
        @test first(h) == p[2]

        p[2].infection_time = 24
        update!(h, 2)
        @test length(h) == 7
        @test first(h).id == 1
        @test first(h) == p[1]

        p[2].infection_time = 8
        push!(h, p[2])
        @test length(h) == 7
        @test first(h).id == 2
        @test first(h) == p[2]

        p[2].infection_time = 24
        push!(h, p[2])
        @test length(h) == 7
        @test first(h).id == 1
        @test first(h) == p[1]

        i = pop!(h)
        @test i.id == 1
        @test i == p[1]
        @test first(h).id == 7

        delete!(h, 3)
        @test !(3 in h)
        @test length(h) == 5

        delete!(h, p[5])
        @test !(5 in h)
        @test length(h) == 4

        @test pop!(h).id == 7
        @test pop!(h).id == 6
        @test pop!(h).id == 2
        @test pop!(h).id == 4
        @test isempty(h)

    end

end

