abstract type StateType{T <: Integer} <: Enum{T} end

function exposed end
function infectious end

"""
    @states StateTypeName[::BaseType] state1[=x] state2[=y]

Create a `StateType{BaseType}` subtype with name `StateTypeName` and state type values of
`value1` and `value2` with optional assigned values of `x` and `y`, respectively.

This follows the same construct as `@enum` with notable differences. It does
not export the members as `const` to avoid polluting the environment. Second,
it an additional method for fetching and converting from `Symbol` to
`StateType{BaseType}`.
"""
macro states(T::Union{Symbol,Expr}, syms...)
# most of the macro is copied from
# https://github.com/JuliaLang/julia/blob/master/base/Enums.jl
    if isempty(syms)
        throw(ArgumentError("no arguments given for StateType $T"))
    end
    basetype = Int32
    typename = T
    if isa(T, Expr) && T.head === :(::) && length(T.args) == 2 && isa(T.args[1], Symbol)
        typename = T.args[1]
        basetype = Core.eval(__module__, T.args[2])
        if !isa(basetype, DataType) || !(basetype <: Integer) || !isbitstype(basetype)
            throw(ArgumentError("invalid base type for StateType $typename, $T=::$basetype; base type must be an integer primitive type"))
        end
    elseif !isa(T, Symbol)
        throw(ArgumentError("invalid type expression for enum $T"))
    end
    values = Vector{basetype}()
    seen = Set{Symbol}()
    namemap = Dict{basetype,Symbol}()
    valuemap = Dict{Symbol,basetype}()
    exposed = Set{basetype}()
    infectious = Set{basetype}()
    lo = hi = 0
    i = zero(basetype)
    hasexpr = false
    isexposed = false
    isinfectious = false

    if length(syms) == 1 && syms[1] isa Expr && syms[1].head === :block
        syms = syms[1].args
    end
    for s in syms
        s isa LineNumberNode && continue
        if isa(s, Symbol)
            if i == typemin(basetype) && !isempty(values)
                throw(ArgumentError("overflow in value \"$s\" of StateType $typename"))
            end
        elseif isa(s, Expr) &&
               (s.head === :(=) || s.head === :kw) &&
               length(s.args) == 2 && isa(s.args[1], Symbol)
            i = Core.eval(__module__, s.args[2]) # allow exprs, e.g. uint128"1"
            if !isa(i, Integer)
                throw(ArgumentError("invalid value for StateType $typename, $s; values must be integers"))
            end
            i = convert(basetype, i)
            s = s.args[1]
            hasexpr = true
        elseif isa(s, Expr) &&
            (s.head === :ref) &&
            length(s.args) > 1 && isa(s.args[1], Symbol)
            for condition in s.args[2:end]
                if condition === :(:exposed)
                    isexposed = true
                elseif condition === :(:infectious)
                    isinfectious = true
                else
                    throw(ArgumentError("invalid condition for StateType $typename, $s; allowed conditions [:exposed, :infectious]"))
                end
            end
            s = s.args[1]
        else
            throw(ArgumentError(string("invalid argument for StateType ", typename, ": ", s)))
        end
        s = s::Symbol
        if !Base.isidentifier(s)
            throw(ArgumentError("invalid name for StateType $typename; \"$s\" is not a valid identifier"))
        end
        if hasexpr && haskey(namemap, i)
            throw(ArgumentError("both $s and $(namemap[i]) have value $i in StateType $typename; values must be unique"))
        end
        namemap[i] = s
        # the valuemap is used by the additional helpers
        valuemap[s] = i
        push!(values, i)
        if s in seen
            throw(ArgumentError("name \"$s\" in StateType $typename is not unique"))
        end
        push!(seen, s)
        if length(values) == 1
            lo = hi = i
        else
            lo = min(lo, i)
            hi = max(hi, i)
        end
        if isexposed
            push!(exposed, i)
        end
        if isinfectious
            push!(infectious, i)
        end
        i += oneunit(i)
        isinfectious = false
        isexposed = false
    end
    blk = quote
        # enum definition
        Base.@__doc__(primitive type $(esc(typename)) <: StateType{$(basetype)} $(sizeof(basetype) * 8) end)
        function $(esc(typename))(x::Integer)
            $(Base.Enums.membershiptest(:x, values)) || Base.Enums.enum_argument_error($(Expr(:quote, typename)), x)
            return bitcast($(esc(typename)), convert($(basetype), x))
        end
        Base.Enums.namemap(::Type{$(esc(typename))}) = $(esc(namemap))
        Base.typemin(x::Type{$(esc(typename))}) = $(esc(typename))($lo)
        Base.typemax(x::Type{$(esc(typename))}) = $(esc(typename))($hi)
        let insts = (Any[ $(esc(typename))(v) for v in $values ]...,)
            Base.instances(::Type{$(esc(typename))}) = insts
        end
        # additional helpers
        function $(esc(typename))(x::Symbol)
            (x in keys($(esc(valuemap)))) || Base.Enums.enum_argument_error($(Expr(:quote, typename)), x)
            return $(esc(typename))($(esc(valuemap))[x])
        end
        Base.convert(::Type{$(esc(typename))}, x::Symbol) = $(esc(typename))(x)
        Tseir.exposed(x::$(esc(typename))) = $(basetype)(x) in $(esc(exposed))
        Tseir.infectious(x::$(esc(typename))) = $(basetype)(x) in $(esc(infectious))
    end
    push!(blk.args, :nothing)
    blk.head = :toplevel
    return blk
end

mutable struct State{T <: StateType}
    type::T
    time::Int32
    location::Int32
    source::Int32
    migration::Int32
end

type(S::State) = S.type

Base.time(S::State) = S.time
Base.time(::Nothing) = typemax(Int32)

location(S::State) = S.location
location(::Nothing) = typemax(Int32)

source(S::State) = S.source
source(::Nothing) = typemax(Int32)

migration(S::State) = S.migration
migration(::Nothing) = typemax(Int32)

function Base.show(io::IO, s::State{T}) where {T <: StateType}
    print(io,
      "State{", T,
      "}(type=", s.type,
      ", time=", s.time,
      ", location=", s.location,
      ", source=", s.source,
      ", migration=", s.migration, ")"
   )
end

function Base.:(==)(x::State, y::State)
    out = (
      (x.type == y.type) &&
      (x.time == y.time) &&
      (x.location == y.location) &&
      (x.source == y.source) &&
      (x.migration == y.migration)
   )
    return out
end

function State(type::T) where {T <: StateType}
    return State(type, typemax(Int32), typemax(Int32), typemax(Int32), typemax(Int32))
end

function State(type::T, time::Int32) where {T <: StateType}
    return State(type, time, typemax(Int32), typemax(Int32), typemax(Int32))
end

function State(type::T, time::Int32, source::Int32) where {T <: StateType}
    return State(type, time, typemax(Int32), source, typemax(Int32))
end

function assign_event_location!(state::State, transitions::Array{Interval})

    t = time(state)

   # binary search on the transition list to find the first
   # coordinate set that the individual visit after getting sick.
    lo = 1; mid = lo; hi = length(transitions) + 1

    while lo < hi
        mid = (lo + hi) >> 1
        if (transitions[mid]._end >= t)
            hi = mid
        else
            lo = mid + 1
        end
    end

    if hi > length(transitions)
        return
    end

    state.location = transitions[hi]._coord
    state.migration = transitions[hi]._end

end
