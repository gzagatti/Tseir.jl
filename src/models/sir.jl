@states StateSIR S I[:exposed, :infectious] R

SIR = Model{StateSIR}()

add_path!(SIR, :S, :I)
add_path!(SIR, :I, :R)
