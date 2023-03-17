using Gmsh: gmsh
using GR 
using LinearAlgebra
using Plots
using LaTeXStrings



#gmsh.fltk.run()


#..1/11 Generate the mesh




#=
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("t1")
lc = 1e-1
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(1., 0,  0, lc, 2)
gmsh.model.geo.addPoint(1., 1., 0, lc, 3)
gmsh.model.geo.addPoint(0, 1., 0, lc, 4)
gmsh.model.geo.addLine(1, 2, 5)
gmsh.model.geo.addLine(2, 3, 6)
gmsh.model.geo.addLine(3, 4, 7)
gmsh.model.geo.addLine(4, 1, 8)
gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 9)
gmsh.model.geo.addPlaneSurface([9], 10)
gmsh.model.setPhysicalName(2, 11, "My surface")
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
if (false) gmsh.write("t1.msh") end 
=#




#..2/11 Get and sort the mesh nodes
#..Observe that although the mesh is two-dimensional,
#..the z-coordinate that is equal to zero is stored as well.
#..Observe that the coordinates are stored contiguously for computational
#..efficiency



#..6/11 initialize global matrix A and global vector f
#..observe that for simplicity we use dense matrix here

# function

function func(x, y)
  v=x+100; w=y-560;
  val = (1000-v*v/30-w*w/10)/10000;
  #val = 1.2^x;
  if val<0
    val=0;
  end
  return 1
end

const mval = 1;

function assembly(meshname,func)
  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)
  gmsh.model.add("t1")
  gmsh.merge(meshname)
  



  node_ids, node_coord, _ = gmsh.model.mesh.getNodes()
  nnodes = length(node_ids)
  #..sort the node coordinates by ID, such that Node one sits at row 1
  tosort = [node_ids node_coord[2:3:end] node_coord[3:3:end]];
  sorted = sortslices(tosort , dims = 1);
  node_ids = sorted[:,1]
  xnode = sorted[:,2]
  ynode = sorted[:,3]

  #..3/11 Plotting the mesh
  if (false) gmsh.fltk.run() end 


  #-----------------------------------



  #..4/11 Get the mesh elements
  #..observe that we get all the two-dimensional triangular elements from the mesh
  element_types, element_ids, element_connectivity = gmsh.model.mesh.getElements(2)
  nelements = length(element_ids[1])
  println(nelements)




  M = zeros(nnodes,nnodes);
  f2 = zeros(nnodes,1)
  K = zeros(nnodes, nnodes)
  #..7/11 Perform a loop over the elements
  for element_id in 1:nelements

    #....retrieve global numbering of the local nodes of the current element
    node1_id = element_connectivity[1][3*(element_id-1)+1]
    node2_id = element_connectivity[1][3*(element_id-1)+2]
    node3_id = element_connectivity[1][3*(element_id-1)+3]

    #....retrieve the x and y coordinates of the local nodes of the current element
    xnode1 = xnode[node1_id]; xnode2 = xnode[node2_id]; xnode3 = xnode[node3_id];
    ynode1 = ynode[node1_id]; ynode2 = ynode[node2_id]; ynode3 = ynode[node3_id];

    #....compute surface area of the current element
    x12 = xnode2 - xnode1; x13 = xnode3-xnode1;
    y12 = ynode2 - ynode1; y13 = ynode3-ynode1;
    area_id = x12*y13 - x13*y12; area_id = abs(area_id)/2

    # compute the v functions on the element
    # a1 * x1 + b1 * y1 + c1 = 1
    # a1 * x2 + b1 * y2 + c1 = 0
    # ...
    # a3 * x3 + b3 * y3 + c3 = 0

    # (x1 0 0 y1 0 0 1 0 0 )
    # (x2 0 0 y2 0 0 1 0 0 )
    # ...
    # (0 0 x3 0 0 y3 0 0 1)

    # unknowns (a1 a2 a3 b1 b2 b3 c1 c2 c3)'
    # rhs (1 0 0 0 1 0 0 0 1)
    lhs = [xnode1 0 0 ynode1 0 0 1 0 0;
      xnode2 0 0 ynode2 0 0 1 0 0;
      xnode3 0 0 ynode3 0 0 1 0 0;
      0 xnode1 0 0 ynode1 0 0 1 0;
      0 xnode2 0 0 ynode2 0 0 1 0;
      0 xnode3 0 0 ynode3 0 0 1 0;
      0 0 xnode1 0 0 ynode1 0 0 1;
      0 0 xnode2 0 0 ynode2 0 0 1;
      0 0 xnode3 0 0 ynode3 0 0 1]
    rhs = [1; 0; 0; 0; 1; 0; 0; 0; 1];
    cof = lhs\rhs;
     ids = [node1_id  node2_id  node3_id]
     cos = [xnode1 ynode1;
      xnode2 ynode2;
      xnode3 ynode3]

     for i = 1:3
      for j = 1:3
        K[ids[i], ids[j]]+= area_id * (cof[0+i]*cof[0+j]+cof[3+i]*cof[3+j])
      end
      f2[ids[i]]+=area_id * (1/3) * func(cos[i,1], cos[i,2]);
      M[ids[i],ids[i]]+=area_id * (1/3);
     end
    
  end
  α = 1;
  C=α*M;
  
  node_ids2, node_coord, _ = gmsh.model.mesh.getNodes(0,2)
  node_idsn = gmsh.model.mesh.getNodesForPhysicalGroup(1, 6)
  node_ids = node_idsn[1]

  bnd_node_ids = union(node_ids,node_ids2)
  K[bnd_node_ids,:] .= 0;
  K[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)))
  f2[bnd_node_ids] .= 0;
  return K,C,M,f2,xnode,ynode
end

#..8/11 Handle the boundary conditions
#..retrieve boundary nodes by loop over corner point and boundary edges


#..9/11 Compute the numerical solution
K,C,M,f2,xnode,ynode = assembly("door_2d.msh", func)
u = K\f2
println("done")

#..10/11 Plot the numerical solution
# tricont(xnode,ynode,u)
trisurf(xnode,ynode,u)
gmsh.finalize()