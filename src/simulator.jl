import DataStructures: BinaryMinHeap
import RandomNumbers: gen_seed

using HDF5
using ProgressBars
using Distributions: Exponential, Geometric
using Random: AbstractRNG

"""
   sir_save(path::AbstractString, results::Dict{Tuple{Int32, Int32, Int32}, Float64},
      transmission_probability, recovery_rate, N::Integer, seed::Integer)

Save the simulation results as an HDF5 file.
"""
function sir_save(filename, results::Dict{Symbol, Dict{K, Float64} where K},
   p::Population, transmission_probability, recovery_rate, N::Integer, seed::Integer)

   infection_time = Array{Int32, 1}(undef, length(results[:infection]))
   infection_source = Array{Int32, 1}(undef, length(results[:infection]))
   infection_location = Array{Int32, 1}(undef, length(results[:infection]))
   infection_amount = Array{Float64, 1}(undef, length(results[:infection]))

   for (ix, (k,v)) in enumerate(results[:infection])
      infection_time[ix], infection_source[ix], infection_location[ix] = k
      infection_amount[ix] = v
   end

   recovery_time = Array{Int32, 1}(undef, length(results[:recovery]))
   recovery_amount = Array{Float64, 1}(undef, length(results[:recovery]))

   for (ix, (k, v)) in enumerate(results[:recovery])
      recovery_time[ix] = k
      recovery_amount[ix] = v
   end

   h5open(filename, "w") do file

      write(file, "infection/time", infection_time)
      write(file, "infection/source", infection_source)
      write(file, "infection/location", infection_location)
      write(file, "infection/amount", infection_amount)

      write(file, "recovery/time", recovery_time)
      write(file, "recovery/amount", recovery_amount)

      write(file, "param/transmission_probability", transmission_probability)
      write(file, "param/recovery_rate", recovery_rate)
      write(file, "param/iter", N)
      write(file, "param/seed", seed)
      write(file, "param/population_size", length(p))

   end

end


"""
   sir_sweep!(p::Population, epidemics_start::Number, epidemics_end::Number,
      transmission_probability::Number, recovery_rate::Number, N::Integer, seed::Integer,
      output_path::AbstractString)

Runs `N` sir simulation shots in population `p`.

For each simulation the disease can strike patient zero any time between
`epidemics_start` and `epidemics_end`.

The random state is initialized according to `seed`.
"""
function sir_sweep!(p::Population, epidemics_start::Number, epidemics_end::Number,
   transmission_probability::Number, recovery_rate::Number, N::Integer, r::AbstractRNG)

   transmission_distribution = Geometric(transmission_probability)
   recovery_distribution = Exponential(1/recovery_rate)

   results = Dict(
      :infection => Dict{Tuple{Int32, Int32, Int32}, Float64}(),
      :recovery => Dict{Int32, Float64}()
   )
   # every 10 minutes
   interval = Int32(600)

   iter = ProgressBar(1:N)
   set_description(iter, "SIR sweep:      ")
   for i in iter
      i0 = rand(r, p)
      t0 = rand(r, Int32(epidemics_start):Int32(epidemics_end))
      sir!(p, i0, t0, transmission_distribution, recovery_distribution, r)
      collect_results!(results, p, t0, interval, 1/N)
   end

   reset!(p)

   return results

end


"""
   sir!(p::Population, i0::Individual, t0::Int32, transmission_distribution::Sampleable,
      recovery_distribution::Sampleable, r::AbstractRNG)

Runs a single sir simulation shot given patient `i0` was infected at time `t0`.
The `transmission_distribution` is the distribution of duration that two
individuals have to be in contact until contagion occurs. The
`recovery_distribution` is distribution of the duration of the disease.
"""
function sir!(p::Population, source::Individual, t0::Int32,
   transmission_distribution::Sampleable, recovery_distribution::Sampleable, r::AbstractRNG)

   infected = PopulationHeap(Base.By(i->i.infection_time))

   # reset all infection parameters in the population
   reset!(p)

   # get and infect the source
   infect!(source, t0, typemax(Int32))
   push!(infected, source)

   # keeping infecting while there are individuals in the heap
   while length(infected) > 0
      print("\r")
      print("Heap size: $(length(infected))")
      infect_others!(infected, p, transmission_distribution, recovery_distribution, r)
   end

end


"""
   infect_others!(h::PopulationHeap, p::Population, transmission_distribution::Sampleable,
      recovery_distribution::Sampleable, r::AbstractRNG)

Determines which contacts individual `i` at the top of the infected heap (ie.
with the earliest infection time out of the non-recovered individuals) will be
able to infect, updating the infected heap accordingly.
"""
function infect_others!(h::PopulationHeap, p::Population,
   transmission_distribution::Sampleable, recovery_distribution::Sampleable, r::AbstractRNG)

   i = pop!(h)
   recover!(i, recovery_distribution, r)

   if i.infection_time == i.recovery_time
      return
   end

   # we forward the model clock up to individual `i` infection time.
   earliest = i.infection_time

   for (otherid, contact_events) in contacts(i)
      j = p[otherid]
      # contact `j` will only have recovered if it has already been popped
      # from the infected heap. In other words, if its recovery time is earlier
      # than `i` infection time
      if j.status != R
         t = transmission(contact_events, earliest, transmission_distribution, r)
         if t == nothing continue end
         # contact `j` can only get infected by `i` before `i` recovers and
         # before it gets infected by anyone else.
         if ((t <= i.recovery_time) & (t <= j.infection_time))
            s = transmission_source(t, i)
            infect!(j, t, s)
            push!(h, j)
         end
      end
   end

end

"""
   transmission(contact_events::Array{ContactEvent}, earliest::Int32,
      transmission_distribution::Sampleable)

Sample the infection contagion time later than a given starting point. The
starting point should be the time that an infected individual transmiting the
disease first got infected.
"""
function transmission(contact_events::Array{ContactEvent}, earliest::Int32,
   transmission_distribution::Sampleable, r::AbstractRNG)
   # binary search on the contact events to find the first event that occurs
   # after the individual getting sick
   lo = 1; mid = lo; hi = length(contact_events) + 1

   while lo < hi
      mid = (lo + hi) >> 1
      if (contact_events[mid].event_end >= earliest)
         hi = mid
      else
         lo = mid + 1
      end
   end

   if hi > length(contact_events)
      return
   end

   # we compute the time the individual got sick in terms of the cumulative
   # time spent together between the individual and the contact.
   # the transmision will then occur in the second following the amount of time
   # spent together.
   transmission_cum_time = (
      contact_events[hi].cum
      - (
         contact_events[hi].event_end
         - max(earliest, contact_events[hi].event_start)
      )
      + floor(Int32, rand(r, transmission_distribution))
      + Int32(1)
   )

   # binary search to find wich contact event led to the infection
   lo = 1; mid = lo; hi = length(contact_events) + 1

   while lo < hi
      mid = (lo + hi) >> 1
      if (contact_events[mid].cum >= transmission_cum_time)
         hi = mid
      else
         lo = mid + 1
      end
   end

   if hi > length(contact_events)
      return
   end

   # convert the strike time from cumulative time spent together to the model's
   # clock time
   transmission_clock_time = contact_events[hi].event_end - (contact_events[hi].cum - transmission_cum_time)

   return transmission_clock_time

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
   if t <= i.migration_time
      return i.infection_source
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
      if (transition_list[mid].event_end >= t)
         hi = mid
      else
         lo = mid + 1
      end
   end

   return transition_list[hi-1].coordset

end

"""
   collect_results!(results::Dict{Tuple{Int32, Int32,Int32},Float64},
      p::Population, t0::Int32, interval::Integer, weight::Integer)

Collect the results of a single sir run on a dictionary of results based on the
final state of population `p`. It records the number of new infections at
regular intervals starting from `t0`. The total number of records is added to
the results with a given `weight`.
"""
function collect_results!(results::Dict{Symbol, Dict{K, Float64} where K},
   p::Population, t0::Int32, interval::Int32, weight::Float64)
   # round to nearest interval from the bottom, eg. if t0 is 16:17:35 and we
   # have an interval equal to 10 minutes then we round t0 to 16:10:00
   t0 = t0 - (t0%interval)
   for i in p
      if i.status == R
         # round to the nearest interval from the top, eg. if infection_time is
         # 16:17:35 and we have an interval equal to 10 minutes then we round
         # the infection_time to 16:20:00
         infection_time = i.infection_time
         infection_time = (infection_time + interval - (infection_time%interval)) - t0

         infection_location = i.infection_location
         infection_source = i.infection_source

         if !haskey(results[:infection], (infection_time, infection_source, infection_location))
            results[:infection][(infection_time, infection_source, infection_location)] = 1 * weight
         else
            results[:infection][(infection_time, infection_source, infection_location)] += 1 * weight
         end

         recovery_time = i.recovery_time
         recovery_time = (recovery_time + interval - (recovery_time%interval)) - t0

         if !haskey(results[:recovery], recovery_time)
            results[:recovery][recovery_time] = 1 * weight
         else
            results[:recovery][recovery_time] += 1 * weight
         end

      end
   end
   return results
end
