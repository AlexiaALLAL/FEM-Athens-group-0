using Gmsh: gmsh
using MeshIO
using FileIO

gmsh.initialize()

mesh = gmsh.model.add("data/VenturaAccelerometer/test_mesh_01.msh")

coords = []
nodes = []
tris = []
