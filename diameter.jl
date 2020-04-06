using LightGraphs
using BALUtils
using StaticGraphs

ba = readbal(ARGS[1])
dia = diameter(StaticGraph(Graph(visibility_graph(ba))))
open(ARGS[2], "w") do io
    write(io, "diameter\n$dia\n")
end
