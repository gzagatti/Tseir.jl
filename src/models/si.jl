@states StateSI S I[:exposed, :infectious]

SI = Model{StateSI}()

add_path!(SI, :S, :I)

function infect!(i::Individual{StateSI}, t::Number, source::Number)
    i.state = State(StateSI(:I), t, source)
    i.infection = i.state
    assign_event_location!(infection(i), transitions(i))
end
