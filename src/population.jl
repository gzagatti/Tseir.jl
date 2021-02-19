
"""
   Population

A population is a set of individuals. We can perform basic set operations
by passing either an instance of an Individual or its id.
"""
struct Population{T} <: AbstractSet{Individual{T<:StateType}}
   dict::Dict{Int32,Individual}

   Population() = new(Dict{Int32,Individual}())
   Population(p::Array{Individual}) = new(Dict{Int32, Individual}(i.id => i for i in p))
   Population(p::Population) = new(Dict{Int32,Individual}(i.id => i for i in p)) 
end

Base.isempty(p::Population) = isempty(p.dict)
Base.length(p::Population) = length(p.dict)

Base.in(i::Individual, p::Population) = haskey(p.dict, i.id)
Base.in(id::Number, p::Population) = haskey(p.dict, Int32(id))

Base.getindex(p::Population, i::Individual) = p.dict[i.id]
Base.getindex(p::Population, id::Number) = p.dict[Int32(id)]

Base.push!(p::Population, i::Individual) = (p.dict[i.id] = i; p)

Base.pop!(p::Population, i::Individual) = pop!(p.dict, i.id)
Base.pop!(p::Population, id::Number) = pop!(p.dict, Int32(id))
Base.pop!(p::Population, i::Individual, default) = (i in p ? pop!(p, i) : default)
Base.pop!(p::Population, id::Number, default) = (id in p ? pop!(p, Int32(id)) : default)


function Base.pop!(pop::Population)
   isempty(pop) && throw(ArgumentError("set must be non-empty"))
   return pop!(pop.dict)[2]
end

Base.delete!(p::Population, i::Individual) = (delete!(p.dict, i.id); p)
Base.delete!(p::Population, id::Number) = (delete!(p.dict, Int32(id)); p)

Base.copy(p::Population) = Base.copymutable(p)

Base.sizehint!(p::Population, s) = (sizehint!(p.dict, s); p)
Base.empty!(p::Population) = (empty!(p.dict); p)
Base.rehash!(p::Population) = (Base.rehash!(p.dict); p)

Base.iterate(p::Population, i...) = iterate(Base.ValueIterator(p.dict), i...)

reset!(p::Population) = for i in p reset!(i) end

"""
   Eager initialize population

Initializes a population set based on two PostgresSQL statements: the first a
list of individual transitions and the second returns a list of contact events.

This is an eager initialization as all contacts and transitions are pre-loaded.
"""
function eager_init_population(transition_stmt::LibPQ.Statement, contact_stmt::LibPQ.Statement)
   p = Population()
   iter = ProgressBar(LibPQ.execute(contact_stmt))
   set_description(iter, "Filling contacts:")
   for contact_details in iter
      id, otherid, coord, event_start, event_end = contact_details
      if !(id in p) 
         i = Individual(id, transition_stmt)
         push!(p, i)
      end
      if !(otherid in p) 
         other = Individual(otherid, transition_stmt)
         push!(p, other)
      end
      push_contact!(p[id], otherid, coord, event_start, event_end)
   end
   iter = ProgressBar(p)
   set_description(iter, "Sorting contacts: ")
   for i in iter
      sort_contacts!(i)
   end
   return p
end
"""
   Lazy initialize population

Initializes a population set based on three PostgresSQL statements: the first
returns a list of unique individuals, the second a list of individual
transitions and the third a list of contact events.

This is a lazy initialization as individuals are only initialized with contact
and transition statements. The contact and transition list are only loaded when
needed.
"""
function lazy_init_population(id_stmt::LibPQ.Statement, transition_stmt::LibPQ.Statement, contact_stmt::LibPQ.Statement)
   p = Population()
   for (id,) in LibPQ.execute(id_stmt)
      i = Individual(id, transition_stmt, contact_stmt)
      push!(p, i)
   end
   return p
end

"""
   Population Heap

A population heap to keep track of individuals according to an arbitrary
ordering.

This data structure was largely inspired by the `MutableBinaryHeap`
implementation of the `DataStructures.jl` package, available at:
<https://github.com/JuliaCollections/DataStructures.jl/blob/master/src/heaps/mutable_binary_heap.jl>
"""
mutable struct PopulationHeap{O <: Base.Ordering} <: AbstractMutableHeap{Individual,Int}
   ordering::O
   nodes::Vector{Individual}
   nodemap::Dict{Int32,Int}

   function PopulationHeap(ordering::Base.Ordering)
      nodes = Vector{Individual}()
      nodemap = Dict{Int32,Int}()
      new{typeof(ordering)}(ordering, nodes, nodemap)
   end

end

function Base.show(io::IO, h::PopulationHeap)
   print(io, "PopulationHeap(")
   nodes = h.nodes
   n = length(nodes)
   if n > 0
      print(io, string(nodes[1]))
      for i = 2:n
         print(io, ", $(nodes[i])")
      end
   end
   print(io, ")")
end

Base.isempty(h::PopulationHeap) = isempty(h.nodes)
Base.length(h::PopulationHeap) = length(h.nodes)

Base.in(i::Individual, h::PopulationHeap) = haskey(h.nodemap, i.id)
Base.in(id::Number, h::PopulationHeap) = haskey(h.nodemap, Int32(id))

Base.getindex(h::PopulationHeap, id::Number) = h.nodes[h.nodemap[id]]
Base.getindex(h::PopulationHeap, i::Individual) = h.nodes[h.nodemap[i.id]]

@inline Base.first(h::PopulationHeap) = h.nodes[1]

"""
   push!(h::PopulationHeap, i::Individual)

Push Individual `i` to population heap `h`. If the heap already contains `i`,
then it updates `i` position.
"""
function Base.push!(h::PopulationHeap, i::Individual)
   nodemap = h.nodemap
   if haskey(nodemap, i.id)
      update!(h, i)
      return i
   end
   ordering = h.ordering
   nodes = h.nodes
   nd_id = length(nodes) + 1
   push!(nodes, i)
   nodemap[i.id] = nd_id
   _heap_bubble_up!(ordering, nodes, nodemap, nd_id)
   return i
end

function Base.sizehint!(h::PopulationHeap, s)
   sizehint!(h.nodes, s)
   sizehint!(h.nodemap, s)
   return h
end

"""
   pop!(h::PopulationHeap)

Remove the individual at the top of the heap from the heap returning it.
"""
Base.pop!(h::PopulationHeap) = _binary_heap_pop!(h.ordering, h.nodes, h.nodemap)

"""
   update!(h::PopulationHeap, i::Individual)

Update the position of individual `i` in population heap `h`.
"""
function update!(h::PopulationHeap, i::Individual)
   nodes = h.nodes
   nodemap = h.nodemap
   ordering = h.ordering
   nd_id = nodemap[i.id]
   if nd_id > 1
      p_id = nd_id >> 1
      p = nodes[p_id]
      if Base.lt(ordering, i, p)
         _heap_bubble_up!(ordering, nodes, nodemap, nd_id)
      else
         _heap_bubble_down!(ordering, nodes, nodemap, nd_id)
      end
   else
      _heap_bubble_down!(ordering, nodes, nodemap, nd_id)
   end
end
update!(h::PopulationHeap, id::Number) = update!(h, h.nodes[h.nodemap[id]])

"""
   delete!(h::PopulationHeap, i::Individual)

Remove individual `i` from population heap `h`.
"""
function Base.delete!(h::PopulationHeap, i::Individual)
   nd_id = h.nodemap[i.id]
   _binary_heap_pop!(h.ordering, h.nodes, h.nodemap, nd_id)
   return h
end
function Base.delete!(h::PopulationHeap, id::Number)
   nd_id = h.nodemap[id]
   _binary_heap_pop!(h.ordering, h.nodes, h.nodemap, nd_id)
   return h
end

function _binary_heap_pop!(ord::Base.Ordering, nodes::Vector{Individual},
   nodemap::Dict{Int32,Int}, nd_id::Int=1)

   # extract node
   rt = nodes[nd_id]
   @inbounds delete!(nodemap, rt.id)

   # if node-to-remove is at end, we can just pop it
   # the same applies to 1-element heaps that are empty after removing the last element
   if nd_id == lastindex(nodes)
      pop!(nodes)
   else
     # move the last node to the position of the node-to-remove
      @inbounds nodes[nd_id] = new_rt = nodes[end]
      pop!(nodes)
      @inbounds nodemap[new_rt.id] = nd_id

      if length(nodes) > 1
         if Base.lt(ord, new_rt, rt)
            _heap_bubble_up!(ord, nodes, nodemap, nd_id)
         else
            _heap_bubble_down!(ord, nodes, nodemap, nd_id)
         end
      end
   end

   return rt

end

function _heap_bubble_up!(ord::Base.Ordering, nodes::Vector{Individual}, 
   nodemap::Dict{Int32,Int}, nd_id::Int)

   @inbounds nd = nodes[nd_id]

   swapped = true
   i = nd_id

   while swapped && i > 1 # nd is not root
      p = i >> 1
      @inbounds nd_p = nodes[p]

      if Base.lt(ord, nd, nd_p)
         # move parent downward
         @inbounds nodes[i] = nd_p
         @inbounds nodemap[nd_p.id] = i
         i = p
      else
         swapped = false
      end
   end

   if i != nd_id
      nodes[i] = nd
      nodemap[nd.id] = i
   end

end

function _heap_bubble_down!(ord::Base.Ordering, nodes::Vector{Individual}, 
   nodemap::Dict{Int32,Int}, nd_id::Int)

   @inbounds nd = nodes[nd_id]

   n = length(nodes)
   last_parent = n >> 1

   swapped = true
   i = nd_id

   while swapped && i <= last_parent
      il = i << 1

      if il < n # contains both left and right children

         ir = il + 1

         # determine the better child
         @inbounds nd_l = nodes[il]
         @inbounds nd_r = nodes[ir]

         if Base.lt(ord, nd_r, nd_l)
            # consider right child
            if Base.lt(ord, nd_r, nd)
               @inbounds nodes[i] = nd_r
               @inbounds nodemap[nd_r.id] = i
               i = ir
            else
               swapped = false
            end
         else
         # consider left child
            if Base.lt(ord, nd_l, nd)
              @inbounds nodes[i] = nd_l
              @inbounds nodemap[nd_l.id] = i
              i = il
            else
              swapped = false
            end
         end

      else # contains only left child

         nd_l = nodes[il]
         if Base.lt(ord, nd_l, nd)
            @inbounds nodes[i] = nd_l
            @inbounds nodemap[nd_l.id] = i
            i = il
         else
            swapped = false
         end

      end

   end

   if i != nd_id
      @inbounds nodes[i] = nd
      @inbounds nodemap[nd.id] = i
   end

end
