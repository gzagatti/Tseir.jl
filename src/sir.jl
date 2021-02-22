@states StateSIR S I[:exposed, :infectious] R

SIR = Model{StateSIR}()

add_path!(SIR, :S, :I)
add_path!(SIR, :I, :R)

function infect!(i::Individual{StateSIR}, m::Model{StateSIR}, t::Number, source::Number)
    i.state = State(:I, t, source)
    i.exposure = i.state
end
