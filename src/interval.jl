"""
   Interval

This structure defines a time interval. Attributes include the start and end
time, the total cumulative duration of prior intervals up to and including this
one and a coordinate indicating location.
"""
struct Interval
    _start::Int32
    _end::Int32
    _cum::Int32
    _coord::Int32
end

