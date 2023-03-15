using Gmsh
using Gridap
using GridapGmsh

model = GmshDiscreteModel("../data/VenturaAccelerometer/test_mesh_v2.msh")
elements = model.grid.cell_node_ids
nodes = model.grid.node_coordinates

Areas = zeros(length(elements))
A = 0

push!(nodes, nodes[1])
print(nodes)

for (i, element) in enumerate(elements)
    A1 = 0
    A2 = 0

    A += abs(A1-A2)/2  
    Areas[i] = abs(A1-A2)/2  

end
println(Areas)
