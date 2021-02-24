@states StateSEIR S E[:exposed] I[:infectious] R

SEIR = Model{StateSEIR}()

add_path!(SEIR, :S, :E)
add_path!(SEIR, :E, :I)
add_path!(SEIR, :I, :R)
