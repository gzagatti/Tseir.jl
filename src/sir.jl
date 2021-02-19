@states StateSIR S I R

SIR = Model{StateSIR}()

set_exposed!(SIR, :I)
set_infectious!(SIR, :I)
set_terminal!(SIR, :R)

add_path!(SIR, :S, :I)
add_path!(SIR, :I, :R)

function infect!(i::Individual{StateSIR}, m::Model{StateSIR}, t::Number, source::Number)
   i.state = State(:I, t, source)
   i.exposure = i.state
end
