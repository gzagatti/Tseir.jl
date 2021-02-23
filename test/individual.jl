if !("$(@__DIR__)/../src"   in LOAD_PATH)
    pushfirst!(LOAD_PATH, "$(@__DIR__)/../src")
end

using Test
using Tseir
using Random
using Distributions:DiscreteUniform

r = Random.GLOBAL_RNG

@states Foo foo bar

@testset "Individual" begin

    individuals = [Individual{Foo}(id) for id in 1:2]

    for ix in 1:length(individuals)
        @test individuals[ix].id == ix
        @test individuals[ix].state == State(Foo(0))
        @test !isdefined(individuals[ix], :transition_list)
        @test !isdefined(individuals[ix], :contact_list)
    end

    transition_list = [
       Tseir.Interval(50, 60, typemax(Int32), 1),
       Tseir.Interval(70, 80, typemax(Int32), 2),
       Tseir.Interval(85, 99, typemax(Int32), 1)
    ]
    Tseir.cat_transition!(individuals[1], transition_list)

    @test length(individuals[1].transition_list) == 3

    transition_list = [
       Tseir.Interval(40, 55, typemax(Int32), 1),
       Tseir.Interval(75, 95, typemax(Int32), 3),
    ]
    Tseir.cat_transition!(individuals[2], transition_list)

    @test length(individuals[2].transition_list) == 2

    for i in individuals
        for ix in 2:length(i.transition_list)
            @test i.transition_list[ix - 1]._end < i.transition_list[ix]._start
            @test i.transition_list[ix - 1]._coord != i.transition_list[ix]._coord
        end
    end

    i1 = individuals[1]
    i2 = individuals[2]

    contact_list = [
        (i1, i2.id, 75, 77),
        (i1, i2.id, 50, 55),
        (i1, i2.id, 78, 80),
        (i2, i1.id, 75, 77),
        (i2, i1.id, 78, 80),
        (i2, i1.id, 50, 55),
    ]
    for (i, otherid, interval_start, interval_end) in contact_list
        Tseir.push_contact!(i, otherid, interval_start, interval_end)
    end

    for i in individuals
        Tseir.sort_contacts!(i)
        for ix in 2:length(i.contact_list)
            @test i.contact_list[ix - 1] <= i.contact_list[ix]
        end
    end

    @test haskey(i1.contact_list, i2.id)
    @test haskey(i2.contact_list, i1.id)

    @test i1.contact_list[i2.id][1]._start == 50
    @test i1.contact_list[i2.id][1]._end == 55
    @test i1.contact_list[i2.id][1]._cum == 5

    @test i1.contact_list[i2.id][2]._start == 75
    @test i1.contact_list[i2.id][2]._end == 77
    @test i1.contact_list[i2.id][2]._cum == 7

    @test i1.contact_list[i2.id][3]._start == 78
    @test i1.contact_list[i2.id][3]._end == 80
    @test i1.contact_list[i2.id][3]._cum == 9

    i1.state = State(typemax(Foo), Int32(10), Int32(5))
    i1.infection = i1.state
    @test infection(i1) == State(typemax(Foo), Int32(10), Int32(5))
    @test time(infection(i1)) == Int32(10)
    @test source(infection(i1)) == Int32(5)
    reset!(i1)
    @test i1.state == State(typemin(Foo))
    @test infection(i1) == nothing
    @test time(infection(i1)) == typemax(Int32)

end
