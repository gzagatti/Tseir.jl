"""
   sweep!(p::Population, m::Model, epidemics_start::Number,
      epidemics_end::Number, N::Integer, seed::Integer)

Runs `N` simulation shots in population `p` using model `m`.

For each simulation the disease can strike patient zero any time between
`epidemics_start` and `epidemics_end`.

The random state is initialized according to `seed`.
"""
function sweep!(p::Population, m::Model, epidemics_start::Number,
   epidemics_end::Number, N::Integer, save_interval::Number)

    results = Dict(
      :infection => Dict{Tuple{Int32,Int32,Int32},Float64}(),
      :recovery => Dict{Int32,Float64}()
    )

    iter = ProgressBar(1:N)
    set_description(iter, "Sweep:      ")
    for i in iter
        i0 = rand(m.r, p)
        t0 = rand(m.r, Int32(epidemics_start):Int32(epidemics_end))
        run!(p, m, i0, t0)
        collect_results!(results, p, t0, Int32(save_interval), 1 / N)
    end

    reset!(p)

    return results

end

"""
   run!(p::Population, m::Model, i0::Individual, t0::Int32)

Runs a single simulation shot on population `p` using model `m` given that
patient `i0` was infected at time `t0`.
"""
function run!(p::Population, m::Model, i0::Individual, t0::Int32)

    infected = PopulationHeap(Base.By(i -> time(infection(i))))

    # reset all infection parameters in the population
    reset!(p)

    # get and infect the source
    infect!(i0, t0, typemax(Int32))
    push!(infected, i0)

    # keeping infecting while there are individuals in the heap
    while length(infected) > 0
        infect_others!(infected, p, m)
    end

end



"""
   infect_others!(h::PopulationHeap, p::Population, m::Model)

Determines which contacts in population `p` the individual at the top of the
infected heap `h` (ie.  with the earliest infection time out of the
non-recovered individuals) will be able to infect according to model `m`,
updating the infected heap accordingly.
"""
function infect_others!(h::PopulationHeap, p::Population, m::Model)

    i = pop!(h)

    while !infectious(type(state(i)))
        advance!(i, m)
    end

    advance!(i, m)
    assign_event_location!(infection(i), transitions(i))

    if time(infection(i)) == time(state(i))
        return
    end

    for (otherid, contact_events) in contacts(i)
        j = p[otherid]
        # contact `j` will only have recovered if it has already been popped
        # from the infected heap. In other words, if its recovery time is earlier
        # than `i` infection time
        if can_infect(j, m)
            if type(state(j)) in m.susceptible
                s = state(j)
            elseif type(previous_state(j)) in m.susceptible
                s = previous_state(j)
            else
                continue
            end
            s.time = time(infection(i))
            t = rand(m, s, contact_events)
            if t == nothing continue end
        # contact `j` can only get infected by `i` before `i` recovers and
        # before it gets infected by anyone else.
            if ((t <= time(state(i))) & (t <= time(infection(j))))
                source = transmission_source(t, i)
                if !(type(state(j)) in m.susceptible) rewind!(j) end
                advance!(j, m, t, source)
                push!(h, j)
            end
        end
    end

end

"""
  transmission_source(t::Int32, i::Individual)

Determine the transmission source, given infection time `t`.

If individual `i` did not migrate to another coordinate set before infecting
others, the contagion source will be the same as individual `i` contagion
source. On the other hand, if individual `i` migrated to another coordinate
set, the contagion source will be its previous visited coordinate set.
"""
function transmission_source(t::Int32, i::Individual)
    if t <= migration(infection(i))
        return source(infection(i))
    end

    transition_list = transitions(i)

   # binary search to determine the previously visited location.  Since there
   # was a succesfull transmission, the binary search must return an index
   # within the transmission_list. Second, since the transmission time must be
   # later than the migration time, it must be the case that individual `i`
   # transitioned at least once, thus we start our search in position 2.
    lo = 2; mid = lo; hi = length(transition_list)

    while lo < hi
        mid = (lo + hi) >> 1
        if (transition_list[mid]._end >= t)
            hi = mid
        else
            lo = mid + 1
        end
    end

    return transition_list[hi - 1]._coord

end
