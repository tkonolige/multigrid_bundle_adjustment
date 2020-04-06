using HDF5

out = ARGS[3]

sm = h5read(ARGS[1], "/eigenvalues")
lg = h5read(ARGS[2], "/eigenvalues")
egs = vcat(sm, lg)
min, max = extrema(egs)
open(out, "w") do io
    write(io,"condition_number,smallest,largest\n$(max/min),$min,$max")
end
