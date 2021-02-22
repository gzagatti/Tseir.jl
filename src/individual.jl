"""
   Individual

The individual is the basic unit of our model. We keep track of its infection
status and details regarding its infection event. We also keep track of all
its contact events and all of its transitions between coordinates sets.
"""
mutable struct Individual{T <: StateType}

    id::Int32

    state::State{T}
    infection::Union{State{T},Nothing}

    transition_stmt::LibPQ.Statement
    contact_stmt::LibPQ.Statement

    transition_list::Array{Interval}
    contact_list::Dict{Int32,Array{Interval}}

    function Individual{T}(id) where {T <: StateType}
        self = new{T}()
        self.id = Int32(id)
        self.state = State(typemin(T))
        self.infection = nothing
        return self
    end

    function Individual{T}(id, transition_list::Array{Interval}) where {T <: StateType}
        self = Individual(id)
        cat_transition!(self, transition_list)
    end

    function Individual{T}(id, transition_stmt::LibPQ.Statement) where {T <: StateType}
        self = Individual(id)
        cat_transition!(self, transition_stmt)
        return self
    end

    function Individual{T}(id, transition_stmt::LibPQ.Statement, contact_stmt::LibPQ.Statement) where {T <: StateType}
        self = Individual(id)
        self.transition_stmt = transition_stmt
        self.contact_stmt = contact_stmt
        return self
    end

end

function Base.show(io::IO, i::Individual{T}) where {T <: StateType}
    print(io, "Individual{", T, "}(", i.id, ", ", state(i), " status, ")
    if infection_time(i) < typemax(Int32)
        print(io, infection_time(i), " infection time, ")
    end
    if isdefined(i, :contact_list)
        print(io, length(i.contact_list), " contact(s) and ")
    else
        print(io, " #undef contact(s) and ")
    end
    if isdefined(i, :transition_list)
        print(io, length(i.transition_list), " transition(s))")
    else
        print(io, " #undef transition(s)) ")
    end
end

"""
   push_contact!(i::Individual, otherid::Number, interval_start::Int32, interval_end::Int32)

Push a new contact event to Individual `i` contact list.

Does not compute cumulative contact time.
"""
function push_contact!(i::Individual, otherid::Number, interval_start::Number, interval_end::Number)
    otherid = Int32(otherid)
    if !isdefined(i, :contact_list)
        i.contact_list = Dict{Int32,Array{Interval}}()
    end
    if !haskey(i.contact_list, otherid)
        i.contact_list[otherid] = Array{Interval,1}()
    end
    push!(i.contact_list[otherid], Interval(interval_start, interval_end, typemax(Int32), typemax(Int32)))
end

"""
   sort_contacts!(i::Individual)

Sort individual `i` contact list such that for each contact the `interval_start`
field is sorted in increasing order, then computes cumulative contact time.
"""
function sort_contacts!(i::Individual)
    if !isdefined(i, :contact_list)
        return
    end
    for intervals in values(i.contact_list)
        sort!(intervals, by=i -> i._start)
        cum = intervals[1]._end - intervals[1]._start
        intervals[1] = Interval(intervals[1]._start, intervals[1]._end, cum, typemax(Int32))
        if length(contacts) > 1
            for ix in 2:length(contacts)
                cum = contacts[ix - 1]._cum + (contacts[ix]._end - contacts[ix]._start)
                contacts[ix] = ContactEvent(contacts[ix]._start, contacts[ix]._end, cum, typemax(Int32))
            end
        end
    end
end

"""
   cat_transition!(i::Individual, transition_list::Array{Interval})

Concatenate a transition list to individual `i` transition list.

Assumes that the transition list is sorted according to `event_start` and
transitions are non-overlapping.
"""
function cat_transition!(i::Individual, transition_list::Array{Interval})
    if !isdefined(i, :transition_list)
        i.transition_list = Array{Interval,1}()
    end
    if length(transition_list) > 1
        for ix in 2:length(transition_list)
            @assert transition_list[ix - 1]._end <= transition_list[ix]._end
            @assert transition_list[ix - 1].coord != transition_list[ix]._coord
        end
    end
    i.transition_list = transition_list
end

function cat_transition!(i::Individual, transition_stmt::LibPQ.Statement)
    if !isdefined(i, :transition_list)
        i.transition_list = Array{Interval,1}()
    end
    for transition in LibPQ.execute(transition_stmt, (i.id,))
        push!(i.transition_list, Interval(transition...))
    end
    sort!(i.transition_list, by=i -> i._start)
end

"""
   contacts(i::Individual)

Get individual `i` contact list.
"""
function contacts(i::Individual)
    if !isdefined(i, :contact_list)
        if isdefined(i, :contact_stmt)
            i.contact_list = Dict{Int32,Array{Interval}}()
            for contact_details in LibPQ.execute(i.contact_stmt, (i.id,))
                otherid, interval_start, interval_end = contact_details
                push_contact!(i, otherid, interval_start, interval_end)
            end
            sort_contacts!(i)
        end
    end
    return i.contact_list
end

"""
   transitions(i::Individual)

Get individual `i` transition list.
"""
function transitions(i::Individual)
    if !isdefined(i, :transition_list)
        if isdefined(i, :transition_stmt)
            cat_transition!(i, i.transition_stmt)
        end
    end
    return i.transition_list
end

"""
   state(i::Individual)

Get individual `i` current state.
"""
function state(i::Individual)
    return i.state
end

"""
   infection(i::Individual)

Get individual `i` exposed state.
"""
function infection(i::Individual)
    return i.infection
end

"""
   can_infect(i::Individual, m::Model)

Returns wether individual `i` can be infected.
"""
function can_infect(i::Individual{T}, m::Model{T}) where {T <: StateType}
    return i.state in m.susceptible
end

"""
   reset!(i::Individual{T}) where {T<:StateType}

Reset individual `i` infection status to the first state and reset all
infection parameters.
"""
function reset!(i::Individual{T}) where {T <: StateType}
    i.state = State(typemin(T))
    i.infection = nothing
end

"""
   infect!(i::Individual{T}, m::Model{T}, t::Number, source::Number)

Infect individual `i`. What an infection is depends on the `StateType` which
means the function needs to be defined for each `StateType`.
"""
function infect!(i::Individual{T}, m::Model{T}, t::Number, source::Number) where {T <: StateType} end

"""
   advance!(i::Individual, m::Model, [t::Number, s::Number])

Advance individual `i` state given model `m`.

If no time `t` and `source`, the next state time is sampled accordingly.
"""
function advance!(i::Individual, m::Model)
    i.state = advance(m, state(i))
    if exposed(type(i.state))
        i.infection = i.state
    end
end

function advance!(i::Individual, m::Model, t::Number, s::Number)
    i.state = advance(m, state(i), t, s)
    if exposed(type(i.state))
        i.infection = i.state
    end
end

