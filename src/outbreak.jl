"""
    Outbreak

Holds information about a disease outbreak. It includes information about
individual states during the outbreak.
"""
struct Outbreak{T <: StateType}
    root_state::State{T}
    states::Dict{Int32,Vector{State{T}}}
    infections::Dict{Int32,Int}

    function Outbreak{T}() where {T <: StateType}
        return new{T}(
            State(typemin(T)),
            Dict{Int32,Vector{State{T}}}(),
            Dict{Int32,Int}()
        )
    end

end

function Base.show(io::IO, o::Outbreak)
    print(io, "Outbreak{", basetype(o), "}(")
    print(length(o.states), " affected and ", length(o.infections), " infected)")
end

Base.isempty(o::Outbreak) = isempty(o.states)
Base.length(o::Outbreak) = length(o.states)

Base.in(i::Individual, o::Outbreak) = haskey(o.states, i.id)
Base.in(id::Number, o::Outbreak) = haskey(o.states, Int32(id))

Base.getindex(o::Outbreak, i::Individual) = i in o ? o.states[i.id] : [copy(o.root_state)]
Base.getindex(o::Outbreak, id::Number) = id in o ? o.states[Int32(id)] : [copy(o.root_state)]

Base.delete!(o::Outbreak, i::Individual) = (delete!(o.states, i.id); delete!(o.infections, i.id); o)
Base.delete!(o::Outbreak, id::Number) = (delete!(o.states, Int32(id)); delete!(o.infections, Int32(id)); o)

Base.copy(o::Outbreak) = Base.copymutable(o)

Base.sizehint!(o::Outbreak, s) = (sizehint!(o.states, s); sizehint!(o.infections, s); o)
Base.empty!(o::Outbreak) = (empty!(o.states); empty!(o.infections); o)
Base.rehash!(o::Outbreak) = (Base.rehash!(o.states); Base.rehash!(o.infections); o)

Base.iterate(o::Outbreak, i...) = iterate(o.states, i...)

"""
   state(o::Outbreak, i::Individual)
   state(o::Outbreak, id::Number)

Get individual `i` current state in outbreak `o`.
"""
state(o::Outbreak, i::Individual) = o[i][end]
state(o::Outbreak, id::Number) = o[id][end]


"""
   previous_state(o::Outbreak, i::Individual)
   previous_state(o::Outbreak, id::Number)

Get individual `i` previous state in outbreak `o`.
"""
previous_state(o::Outbreak, i::Individual) = (s = o[i]; length(s) == 1 ? copy(o.root_state) : s[length(s) - 1])
previous_state(o::Outbreak, id::Number) = (s = o[id]; length(s) == 1 ? copy(o.root_state) : s[length(s) - 1])

"""
   infection(i::Individual)
   infection(id::Number)

Returns individual `i` infected state in outbreak `o`.
"""
infection(o::Outbreak, i::Individual) = infection(o, i.id)
        function infection(o::Outbreak, id::Number)
    if haskey(o.infections, Int32(id))
        return o.states[Int32(id)][o.infections[Int32(id)]]
    else
        return nothing
    end
end

"""
   susceptible(o::Outbreak, i::Individual, m::Model)
   susceptible(o::Outbreak, id::Number, m::Model)

Returns individual `i` susceptible state in outbreak `o` given model `m`.
"""
susceptible(o::Outbreak, i::Individual, m::Model) = susceptible(o, i.id, m)
function susceptible(o::Outbreak, id::Number, m::Model)
    if type(state(o, id)) in m.susceptible
        return state(o, id)
    elseif type(previous_state(o, id)) in m.susceptible
        return previous_state(o, id)
    else
        return nothing
    end
end

"""
   infect!(o::Outbreak{T}, i::Individual{T}, t::Number, source::Number) where {T <: StateType}

Infect individual `i` in outbreak `o`. What an infection is depends on
`StateType` which means the function needs to be defined for each `StateType`.
"""
function infect!(o::Outbreak{T}, i::Individual, t::Number, source::Number) where {T <: StateType} end

"""
   advance!(o::Outbreak, i::Individual, m::Model, [t::Number, s::Number])
   advance!(o::Outbreak, id::Number, m::Model, [t::Number, s::Number])

Advance individual `i` state in outbreak `o` given model `m`.

If no time `t` and `source`, the next state time is sampled accordingly.
"""
advance!(o::Outbreak, i::Individual, m::Model) = advance!(o, i.id, m)
function advance!(o::Outbreak, id::Number, m::Model)
    new_state = advance(m, state(o, id))
    if !(id in o)
        o.states[id] = [new_state]
    else
        push!(o.states[id], new_state)
    end
    if exposed(type(new_state))
        o.infections[id] = length(o.states[id])
    end
end

advance!(o::Outbreak, i::Individual, m::Model, t::Number, s::Number) = advance!(o, i.id, m, t, s)
function advance!(o::Outbreak, id::Number, m::Model, t::Number, s::Number)
    new_state = advance(m, state(o, id), t, s)
    if !(id in o)
        o.states[id] = [new_state]
    else
        push!(o.states[id], new_state)
    end
    if exposed(type(new_state))
        o.infections[id] = length(o.states[id])
end
end

"""
    rewind!(o::Outbreak, i::Individual)
    rewind!(o::Outbreak, id::Number)

Rewind individual `i` state in outbreak `o`.
"""
rewind!(o::Outbreak, i::Individual) = rewind!(o, i.id)
function rewind!(o::Outbreak, id::Number)
    if (id in o)
        if haskey(o.infections, id) && (o.infections[id] == length(o.states[id]))
            delete!(o.infections, id)
        end
        if length(o.states[id]) > 1
            pop!(o.states[id])
        else
            delete!(o.states, id)
        end
    end
    return
end
