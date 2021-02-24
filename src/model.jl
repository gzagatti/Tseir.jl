"""
   Model{T<:StateType}([r::AbstractRNG])

   Constructs a `Model` of `StateType` with optional random seed `r` for
   replication purposes.

   A model holds information about which states are terminal and infectious,
   the intervent distribution for each state (ie. the distribution of time
   until the next event occur given an event type) and the transition rules
   between state types. These information should be added after model
   initialization. 
"""
mutable struct Model{T <: StateType}
    r::AbstractRNG
    susceptible::Set{T}
    interevent_distribution::Dict{T,Sampleable}
    dag::Dict{T,NamedTuple{(:states, :weights),Tuple{Array{T,1},Weights}}}
end

function Model{T}(r::AbstractRNG) where {T <: StateType}
    susceptible = Set{T}()
    interevent_distribution = Dict{T,Sampleable}()
    dag = Dict{T,NamedTuple{(:states, :weights),Tuple{Array{T,1},Weights}}}()
    return Model(r, susceptible, interevent_distribution, dag)
end

function Model{T}() where {T <: StateType}
    r = Random.default_rng()
    return Model{T}(r)
end

function Base.show(io::IO, m::Model{T}) where {T <: StateType}
    print(io, "Model{", T, "}")
end

basetype(::Model{T}) where {T <: StateType} = T

"""
   set_rng!(m::Model, r::AbstractRNG)

   Set the randon number generator seed for the model.
"""
function set_rng!(m::Model, r::AbstractRNG)
    m.r = r
end

"""
   set_interevent_distribution!(m::Model{T}, s::Union{T, Symbol}, d::Sampleable)

Sets the interevent_distribution of state type `s` to distribution `d` in model `m`.
"""
function set_interevent_distribution!(m::Model{T}, s::Union{T,Symbol}, d::Sampleable) where {T <: StateType}
    s = convert(T, s)
    m.interevent_distribution[s] = d
end


"""
   add_path!(m::Model{T}, source::T, sink::Array{Union{T, Symbol}, 1}, w<:AbstractArray{<:Number, 1}) where {T<:StateType}
   add_path!(m::Model{T}, source::Union{T, Symbol}, sink::Union{T, Symbol}, [w::I]) where {T<:StateType, I<:Number}

Adds the rule to move from `source` to `sink` in model `m`. Weights `w`
determine the relative odds of moving to each state listed in `sink`.
"""
function add_path!(m::Model{T}, source::Union{T,Symbol}, sink::Union{Array{T,1},Array{Symbol,1}}, w::AbstractArray{<:Number,1}) where {T <: StateType}
    source = convert(T, source)
    sink = convert(Array{T,1}, sink)
    w = Weights(w)
    size(sink) == size(w) ? nothing : throw(ArgumentError("sink and weights are not the same size."))
    m.dag[source] = NamedTuple{(:states, :weights),Tuple{Array{T,1},Weights}}((sink, w))
    for s in sink
        if exposed(s)
push!(m.susceptible, source)
            break
        end
    end
end

function add_path!(m::Model{T}, source::Union{T,Symbol}, sink::Union{T,Symbol}, w::I) where {T <: StateType,I <: Number}
    add_path(m, source, [sink], [w])
end

function add_path!(m::Model{T}, source::Union{T,Symbol}, sink::Union{T,Symbol}) where {T <: StateType}
    add_path!(m, source, [sink], [1])
end

"""
   interevent_distribution(m::Model{T}, state::State{T}) where {T<:StateType}

Return the intervent distribution given model `m` and state `s`.
"""
function interevent_distribution(m::Model{T}, s::T) where {T <: StateType}
    return m.interevent_distribution[s]
end

"""
   next_state(m::Model{T}, s::State{T}) where {T<:StateType}

Returns the next type given model `m` rules for moving from state `s`.
"""
function next_state_type(m::Model{T}, s::T) where {T <: StateType}
    (s, weights) = m.dag[s]
    if size(weights) == (1,) 
        return s[1]
    end
    return sample(s, weights)
end

"""
   advance(m::Model{T}, s::State{T}, [t::Int32, source::Int32]) where {T<:StateType}

Returns the next state given model `m` rules for moving from state `s`.

If no time `t` and `source`, the next state time is sampled accordingly.
"""
function advance(m::Model, s::State)
    return State(next_state_type(m, type(s)), rand(m, s))
end

function advance(m::Model, s::State, t::Int32, source::Int32)
    return State(next_state_type(m, type(s)), t, source)
end

struct SamplerModel{T <: StateType} <: Random.Sampler{Float64}
    m::Model{T}
    s::State{T}
    support::Union{Array{Interval,1},Nothing}
end

function Sampler(::Type{<:AbstractRNG}, m::Model{T}, s::State{T}, support::Union{Array{Interval,1},Nothing}, ::Random.Repetition) where {T <: StateType}
    return SamplerModel{T}(m, s, support)
end
Sampler(rng::AbstractRNG, m::Model{T}, s::State{T}, support::Union{Array{Interval,1},Nothing}, r::Random.Repetition) where {T <: StateType} = Sampler(typeof(rng), m, s, support, r)


"""
   rand([rng::AbstractRNG], m::Model, s::State, [support::Array{Interval, 1}])

Samples the next state time from model `m` in state `s` from the `support`.

If no `support`, the next state time is sampled from continuous time.
"""
rand

Base.rand(m::Model{T}, s::State{T}) where {T <: StateType} = rand(m.r, m, s, nothing)
Base.rand(rng::AbstractRNG, m::Model{T}, s::State{T}) where {T <: StateType} = rand(rng, m, s, nothing)
Base.rand(m::Model{T}, s::State{T}, support::Array{Interval,1}) where {T <: StateType} = rand(m.r, m, s, support)
function Base.rand(rng::AbstractRNG, m::Model{T}, s::State{T}, support::Union{Array{Interval,1},Nothing}) where {T <: StateType}
    return rand(rng, Sampler(rng, m, s, support, Val(1)))
end

function Base.rand(rng::AbstractRNG, sp::SamplerModel{<:StateType})

    m = sp.m
    s = sp.s
    support = sp.support
    d = interevent_distribution(m, type(s))

   # if there is no support, the next state time is simply sampled from
   # continuous clock time.
    if isnothing(support)
        return time(s) + floor(Int32, rand(rng, d))
    end

   # binary search on the support to find the first interval that occurs
   # after the state
    lo = 1; mid = lo; hi = length(support) + 1

    while lo < hi
        mid = (lo + hi) >> 1
        if (support[mid]._end >= time(s))
    hi = mid
        else
            lo = mid + 1
        end
    end

    if hi > length(support)
        return
    end

   # we find the state time in terms of the support. The state will occur in
   # the second following the required amount of time elapses in the support.
    state_support_time = (
      support[hi]._cum
      - (
         support[hi]._end
         - max(time(s), support[hi]._start)
      )
      + floor(Int32, rand(m.r, d))
   )

   # binary search to find wich interval led to the infection
    lo = 1; mid = lo; hi = length(support) + 1

    while lo < hi
        mid = (lo + hi) >> 1
        if (support[mid]._cum >= state_support_time)
    hi = mid
        else
            lo = mid + 1
        end
    end

    if hi > length(support)
        return
    end

   # convert the state time in terms of the support to the model's clock time
    state_clock_time = support[hi]._end - (support[hi]._cum - state_support_time)

    return state_clock_time

end

