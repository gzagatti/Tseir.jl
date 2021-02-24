using Test
using Tseir
using Random
using Distributions:DiscreteUniform

r = Random.GLOBAL_RNG

@states Foo foo bar[:exposed]

@testset "Outbreak" begin

    individuals = [Individual(id) for id in 1:2]
    i1 = individuals[1]
    o = Outbreak{Foo}()

    @test isempty(o)
    @test length(o) == 0
    @test o[i1] == [State(typemin(Foo))]
    @test o[1] == [State(typemin(Foo))]
    @test state(o, i1) == State(typemin(Foo))
    @test previous_state(o, i1) == State(typemin(Foo))
    @test isnothing(infection(o, i1))

    infect!(o, i1, Int32(10), Int32(5))

    @test length(o) == 1
    @test i1 in o
    @test 1 in o
    @test o[i1] == [State(typemax(Foo), Int32(10), Int32(5))]
    @test state(o, i1) == State(typemax(Foo), Int32(10), Int32(5))
    @test previous_state(o, i1) == State(typemin(Foo))
    @test infection(o, i1) == State(typemax(Foo), Int32(10), Int32(5))
    @test time(infection(o, i1)) == Int32(10)
    @test source(infection(o, i1)) == Int32(5)

    rewind!(o, i1)
    @test length(o) == 0

end
