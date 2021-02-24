@states StateSI S I[:exposed, :infectious]

SI = Model{StateSI}()

add_path!(SI, :S, :I)
