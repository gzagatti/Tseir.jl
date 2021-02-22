@states StateSIR S I[:exposed, :infectious] R

SIR = Model{StateSIR}()

add_path!(SIR, :S, :I)
add_path!(SIR, :I, :R)

function infect!(i::Individual{StateSIR}, t::Number, source::Number)
    i.state = State(StateSIR(:I), t, source)
    i.infection = i.state
    assign_event_location!(infection(i), transitions(i))
end
