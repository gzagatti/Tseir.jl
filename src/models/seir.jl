@states StateSEIR S E[:exposed] I[:infectious] R

SEIR = Model{StateSEIR}()

add_path!(SEIR, :S, :E)
add_path!(SEIR, :E, :I)
add_path!(SEIR, :I, :R)

function infect!(i::Individual{StateSEIR}, t::Number, source::Number)
    i.state = State(StateSEIR(:E), t, source)
    i.infection = i.state
    assign_event_location!(infection(i), transitions(i))
end
