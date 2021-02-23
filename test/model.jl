if !("$(@__DIR__)/../src"   in LOAD_PATH)
    pushfirst!(LOAD_PATH, "$(@__DIR__)/../src")
end

using Test
using Tseir
using Random
using StatsBase: Weights
using DataStructures:extract_all!
using Distributions:DiscreteUniform

r = Random.GLOBAL_RNG

@states Foo foo bar[:exposed, :infectious]

@testset "Model" begin

    exposed(Foo(0))

    m = Model{Foo}()
    @test m.r == Random.default_rng()

    set_rng!(m, r)
    @test m.r == r

    m = Model{Foo}(r)
    @test m.r == r

    add_path!(m, :foo, :bar)
    @test Foo(:foo) in keys(m.dag)
    @test m.dag[Foo(:foo)] == NamedTuple{(:states, :weights),Tuple{Array{Foo,1},Weights}}(([Foo(:bar)], Weights([1])))
    @test Foo(:foo) in m.susceptible

    set_interevent_distribution!(m, Foo(:foo), DiscreteUniform(1, 1))
    @test m.interevent_distribution[Foo(:foo)] == DiscreteUniform(1, 1)
    @test interevent_distribution(m, Foo(:foo)) == DiscreteUniform(1, 1)

    t = rand(m, State(Foo(:foo), Int32(0)))
    @test t == 1

    support = [
        Tseir.Interval(10, 15, 5, typemax(Int32)),
        Tseir.Interval(17, 20, 8, typemax(Int32)),
        Tseir.Interval(32, 41, 17, typemax(Int32)),
        Tseir.Interval(45, 50, 22, typemax(Int32))
    ]

    set_interevent_distribution!(m, Foo(:foo), DiscreteUniform(1, 1))
    t = rand(r, m, State(Foo(:foo), Int32(5)), support)
    @test t == 11

    set_interevent_distribution!(m, Foo(:foo), DiscreteUniform(2, 2))
    t = rand(r, m, State(Foo(:foo), Int32(18)), support)
    @test t == 20

    set_interevent_distribution!(m, Foo(:foo), DiscreteUniform(11, 11))
    t = rand(r, m, State(Foo(:foo), Int32(19)), support)
    @test t == 46

    set_interevent_distribution!(m, Foo(:foo), DiscreteUniform(2, 2))
    t = rand(r, m, State(Foo(:foo), Int32(15)), support)
    @test t == 19

    set_interevent_distribution!(m, Foo(:foo), DiscreteUniform(1, 1))
    t = rand(r, m, State(Foo(:foo), Int32(52)), support)
    @test t == nothing

    set_interevent_distribution!(m, Foo(:foo), DiscreteUniform(50, 50))
    t = rand(r, m, State(Foo(:foo), Int32(30)), support)
    @test t == nothing

end
