if !("$(@__DIR__)/../src"   in LOAD_PATH)
    pushfirst!(LOAD_PATH, "$(@__DIR__)/../src")
end

using Test
using Tseir
using Random
using Distributions:DiscreteUniform

r = Random.GLOBAL_RNG

@states Foo foo bar

@testset "Population tests" begin 

    @testset "Population" begin

        p = Population{Foo}()

        @test isempty(p) == true
        @test length(p) == 0

        for i in 1:7
            push!(p, Individual{Foo}(i))
        end

        @test isempty(p) == false
        @test length(p) == 7

        @test 1 in p
        @test Individual{Foo}(1) in p

        i = p[2]
        @test typeof(i) == Individual{Foo}
        @test i.id == 2

        i = p[Individual{Foo}(2)]
        @test typeof(i) == Individual{Foo}
        @test i.id == 2

        push!(p, Individual{Foo}(1))
        @test length(p) == 7

        i = pop!(p, Individual{Foo}(1))
        @test i.id == 1
        @test typeof(i) == Individual{Foo}
        @test length(p) == 6

        i = pop!(p, 2)
        @test i.id == 2
        @test typeof(i) == Individual{Foo}
        @test length(p) == 5

        q = delete!(p, Individual{Foo}(3))
        @test q == p
        @test length(p) == 4

        delete!(p, Individual{Foo}(3))
        @test length(p) == 4

        q = delete!(p, 4)
        @test q == p
        @test length(p) == 3

        delete!(p, 4)
        @test length(p) == 3

        i = pop!(p)
        @test typeof(i) == Individual{Foo}
        @test length(p) == 2

        for i in p
            i.state = State(typemax(Foo))
        end

        reset!(p)
        for i in p
            @test i.state == State(typemin(Foo))
        end

        empty!(p)
        @test isempty(p) == true

    end

    @testset "PopulationHeap" begin

        h = PopulationHeap(Base.By(i -> time(infection(i))))
        p = Population{Foo}()

        @test isempty(h) == true
        @test length(h) == 0

        for id in 1:7
            i = Individual{Foo}(id)
            infection_time = floor(Int32, 100 / id)
            i.infection = State(typemax(Foo), infection_time)
            push!(p, i)
            push!(h, i)
        end

        @test isempty(h) == false
        @test length(h) == 7

        @test Individual{Foo}(3) in h
        @test 3 in h

        i = h[Individual{Foo}(2)]
        @test typeof(i) == Individual{Foo}
        @test i.id == 2

        i = h[2]
        @test typeof(i) == Individual{Foo}
        @test i.id == 2
        @test i == p[2]

        i = first(h)
        @test typeof(i) == Individual{Foo}
        @test i.id == 7
        @test i == p[7]

        infection(p[1]).time = 10
        update!(h, p[1])
        @test length(h) == 7
        @test first(h).id == 1
        @test first(h) == p[1]

        infection(p[2]).time = 8
        update!(h, 2)
        @test length(h) == 7
        @test first(h).id == 2
        @test first(h) == p[2]

        infection(p[2]).time = 24
        update!(h, 2)
        @test length(h) == 7
        @test first(h).id == 1
        @test first(h) == p[1]

        infection(p[2]).time = 8
        push!(h, p[2])
        @test length(h) == 7
        @test first(h).id == 2
        @test first(h) == p[2]

        infection(p[2]).time = 24
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
