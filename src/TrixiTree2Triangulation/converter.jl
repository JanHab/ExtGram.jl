#=
Author: Matthias Geratz | Matthias.Geratz@rwth-aachen.de
March 2023


Triangulation of a 2D TreeMesh including solution transfer to the new unique vertices and subsequent output as file.

For assumptions/requirements of the code and for notes on performance see the bottom of this file.
Note: triangle normal definitions might be weird/inconsistent due to the order of the vertices


Useful Trixi files:
https://github.com/trixi-framework/Trixi.jl/blob/main/src/meshes/abstract_tree.jl
https://github.com/trixi-framework/Trixi.jl/blob/main/src/meshes/mesh_io.jl

not used, but seems relevant:
https://github.com/trixi-framework/Trixi.jl/blob/24a03e360ce730aa0fc8b6ddbe509f1a64ef7351/utils/trixi2txt.jl
there is also a to vtk method somewhere.
=#

using Trixi
using LinearAlgebra





# utility function used in generateActiveCellsToVerticesMap()
function is_periodic_neighbor(tree, cell_id, neighbor_id, direction)
  # depends on direction...
  dim = 1 + (direction > 2) # =Int(ceil(direction / 2)) # dimension to check (x or y)
  # if e.g. left neighbor has larger x coordinate than the current cell, its a periodic boundary neighbor
  periodic = (tree.coordinates[dim, neighbor_id] > tree.coordinates[dim, cell_id])
  # invert the result if positive direction (periodic then means not periodic)
  return xor(periodic, iseven(direction))
end


#=
Generates a matrix that stores in every row the unique vertex indices of a given cell.
The cell ordering corresponds to the active cell index of the tree mesh, while the unique
vertex indices are alloted here and do not correspond to other vertex indices/ids.
The ordering is [bottom left, bottom right, top right, top left,
                 left midpoint, right midpoint, bottom midpoint, top midpoint]
An index of zero means the vertex does not exist (possible for the edge midpoints).
Returns active_cells_to_vertices, n_vertices
=#
function generateActiveCellsToVerticesMap(tree)
  active_cell_ids = Trixi.leaf_cells(tree)
  n_active_cells = Trixi.count_leaf_cells(tree)
  active_cells_to_vertices = zeros(Int, n_active_cells, 8)
  next_vertex_index = 1
  # for a given direction, what vertices are the corners of that edge?
  corner_indices = [[1, 4], [2, 3], [1, 2], [4, 3]]
  opposite_direction = [2, 1, 4, 3]

  # generate the inverse mapping: cell_id -> index in active_cell_ids
  cell_id_to_active_idx = zeros(Int, length(tree)) # n_cells = length(mesh.tree)
  for i=1:n_active_cells
    cell_id_to_active_idx[active_cell_ids[i]] = i
  end


  # for all cells: retrieve vertex indices from neighbors or create them if they do not exist yet
  for i=1:n_active_cells
    cell_id = active_cell_ids[i]
    #=
    Only coarse and same resolution neighbors are processed.
    Therefore, to establish coupling with the fine neighbors, the interaction is handled from fine to coarse.
    Information about the coarse neighbors is stored to send newly created vertices to that neighbor if needed.
    In the tree mesh, a cell can have at most 2 coarse neighbors (the other two are same resolution siblings).
    For these neighbors, store: id, position of midpoint (local and neighbor), position of corner (local and neighbor) =#
    coarse_neighbor_info = zeros(Int, 5,2)

    for direction=1:4
      # same resolution neighbor
      neighbor_id = tree.neighbor_ids[direction, cell_id]
      if neighbor_id > 0
        if (Trixi.is_leaf(tree, neighbor_id)
            && neighbor_id < cell_id # neighbor already processed (has corner indices != 0)
            # note: # <=> cell_id_to_active_idx[neighbor_id] == i (active ids monotonically ascending)
            && !is_periodic_neighbor(tree, cell_id, neighbor_id, direction)) # ignore periodic neighbors
          # retrieve indices of the two shared corner vertices
          active_cells_to_vertices[i, corner_indices[direction]] =
          active_cells_to_vertices[cell_id_to_active_idx[neighbor_id], corner_indices[opposite_direction[direction]]]
        end

      # coarse neighbor
      elseif ((neighbor_id = tree.neighbor_ids[direction, tree.parent_ids[cell_id]]) > 0
              && Trixi.is_leaf(tree, neighbor_id)
              && !is_periodic_neighbor(tree, cell_id, neighbor_id, direction))
        # determine which of the two edge corners is corner and which is midpoint of coarser cell edge
        corner_local_pos = 0; # position / index in current cell
        corner_neighbor_pos = 0; # position / index in neighbor cell
        midpoint_local_pos = 0;
        dim = 2 - (direction > 2)
        if tree.coordinates[dim, cell_id] > tree.coordinates[dim, neighbor_id]
          midpoint_local_pos, corner_local_pos = corner_indices[direction]
          corner_neighbor_pos = corner_indices[opposite_direction[direction]][2];
        else
          corner_local_pos, midpoint_local_pos = corner_indices[direction]
          corner_neighbor_pos = corner_indices[opposite_direction[direction]][1];
        end

        # A cell does not retrieve corner indices from fine neighbors,
        # therefore, a cell needs to send its corners to its coarse neighbors.
        
        # store the midpoint to send the index later
        midpoint_neighbor_pos =  4 + opposite_direction[direction];
        neighbor_active_idx = cell_id_to_active_idx[neighbor_id]
        slot = 2 - Int(coarse_neighbor_info[1, 1] == 0) # =1 if first slop is free (entry==0) else = 2
        coarse_neighbor_info[1:3, slot] = [neighbor_active_idx, midpoint_local_pos, midpoint_neighbor_pos]

        # retrieve neighbor corner vertex if it exists, else mark it to send the index later
        neighbor_corner_idx = active_cells_to_vertices[neighbor_active_idx, corner_neighbor_pos];
        if neighbor_corner_idx > 0 # neighbor already has corner idx
          active_cells_to_vertices[i, corner_local_pos] = neighbor_corner_idx # retrieve
        else
          coarse_neighbor_info[4:5, slot] = [corner_local_pos, corner_neighbor_pos]
        end

      # else: neighbor is finer or does not exist (boundary)
      end # if elsif neighbors
    end # for direction


    # assign new indices to the uninitialized corners
    for j=1:4
      if active_cells_to_vertices[i, j] == 0
        active_cells_to_vertices[i, j] = next_vertex_index
        next_vertex_index += 1
      end
    end

    # send indices to the stored coarse neighbors (if needed)
    for j=1:2
      neighbor_active_idx = coarse_neighbor_info[1, j]
      if neighbor_active_idx > 0
        # send midpoint
        local_pos = coarse_neighbor_info[2, j]
        neighbor_pos = coarse_neighbor_info[3, j]
        active_cells_to_vertices[neighbor_active_idx, neighbor_pos] = active_cells_to_vertices[i, local_pos]

        # send corner
        local_pos = coarse_neighbor_info[4, j]
        neighbor_pos = coarse_neighbor_info[5, j]
        if local_pos > 0
          active_cells_to_vertices[neighbor_active_idx, neighbor_pos] = active_cells_to_vertices[i, local_pos]
        end
      end
    end
  end # end for all active cells

  return active_cells_to_vertices, next_vertex_index-1
end







#=
Inverts the mapping i_active_cell->vertices to obtain the mapping of a given i_vertex
to all the cells it is contained within.
The mapping vertex->cells is a matrix where [i_vertex, :] = [9]
Where [1:4] are the cell indices (in active cells) of conected cells
(Vertex maximally part of four cells, if less, then the rest of the four entires is =0).
[5:8] Give the local position (1-8) inside the respective cell ([5] gives pos in [1]).
[9=end] gives the number of cells the vertex is contained within (usefull for iteration and solution averaging).
=#
function generateActiveVerticesToCellsMap(active_cells_to_vertices, n_vertices)
  active_vertices_to_cells = zeros(Int, n_vertices, 9)
  n_active_cells = size(active_cells_to_vertices)[1]

  for i=1:n_active_cells
    for j=1:8
      vertex = active_cells_to_vertices[i, j]
      if vertex > 0
        active_vertices_to_cells[vertex, end] += 1
        active_vertices_to_cells[vertex, active_vertices_to_cells[vertex, end]] = i
        active_vertices_to_cells[vertex, 4 + active_vertices_to_cells[vertex, end]] = j
      end
    end
  end # for all active cells

  return active_vertices_to_cells
end



function collectCoordinates(tree, active_vertices_to_cells)
  n_vertices = size(active_vertices_to_cells)[1]
  coordinates = zeros(n_vertices, 2)

  active_cell_ids = Trixi.leaf_cells(tree)
  # [vertex position] -> vector from center to that vertex on reference cell
  position_lookup = [[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, 0], [1, 0], [0, -1], [0, 1]];

  for i=1:n_vertices
    # get the first cell the vertex is part of an calculate coordinate via that cell
    active_cell_id = active_vertices_to_cells[i, 1]
    position = active_vertices_to_cells[i, 5]

    # get information from Trixi
    global_cell_id = active_cell_ids[active_cell_id]
    cell_length = Trixi.length_at_cell(tree, global_cell_id)
    cell_coordinates = Trixi.cell_coordinates(tree, global_cell_id)
    delta = 0.5*cell_length

    coordinates[i, :] = cell_coordinates + delta .* position_lookup[position]
  end # for all vertices

  return coordinates
end








function triangulateMesh(active_cells_to_vertices)
  n_cells = size(active_cells_to_vertices)[1]

  # the number of triangles in a given cell is 2 + the number of midpoints (each midpoint add a triangle)
  n_triangles_total = 0
  n_triangles_per_cell = zeros(Int, n_cells)
  for i=1:n_cells
    n_triangles_per_cell[i] = sum(active_cells_to_vertices[i, :] .> 0) - 2
    n_triangles_total += n_triangles_per_cell[i]
  end
  connectivity = zeros(Int, 3, n_triangles_total)


  # maps the presence of midpoints [i,j,k,l] to the correct index in the triangulation table
  table_index = reshape([1,2,3,10,4,7,8,12,5,6,9,14,11,13,15,16], 2, 2, 2, 2);
  # for a given table_index gives the indices for the vertices array to yield the connectivity entries
  triangulation_table = [[1 2 3; 1 3 4],
                         [1 2 5; 2 3 5; 3 4 5],
                         [1 2 6; 1 6 4; 3 4 6],
                         [1 7 4; 3 4 7; 2 3 7],
                         [1 8 4; 1 2 8; 2 3 8],
                         [1 2 5; 2 8 5; 2 3 8; 4 5 8],
                         [1 7 5; 2 3 7; 3 5 7; 3 4 5],
                         [1 7 4; 7 6 4; 2 6 7; 3 4 6],
                         [1 2 6; 1 6 8; 1 8 4; 3 8 6],
                         [1 2 6; 1 6 5; 5 6 3; 3 4 5],
                         [1 7 8; 1 8 4; 2 3 7; 3 8 7],
                         [1 7 5; 2 6 7; 7 6 5; 3 5 6; 3 4 5],
                         [1 7 5; 7 8 5; 4 5 8; 2 3 7; 3 8 7],
                         [1 2 6; 1 6 5; 5 6 8; 5 8 4; 3 8 6],
                         [1 7 8; 1 8 4; 7 6 8; 2 6 7; 3 8 6],
                         [1 7 5; 2 6 7; 7 6 5; 5 6 8; 6 3 8; 5 8 4]]


  i_conn = 1
  for i=1:n_cells
    vertices = active_cells_to_vertices[i, :]
    n_triangles = n_triangles_per_cell[i]
    # depending on the vertex signature ([i, :] > 0) of the cell,
    # select the triangluation of the cell, substitute the vertice position with their index and add to connectivity matrix
    connectivity[:, i_conn:i_conn+n_triangles-1] = vertices[triangulation_table[table_index[((vertices[5:8] .> 0) .+1)...]]]'
    i_conn += n_triangles
  end

  return connectivity
end





# evaluates the Lagrange basis function for a set of roots at a point x
function lagrangeBasisEvaluation(n_roots, roots, x)
  l = ones(n_roots)
  for j=1:n_roots
    for i=1:n_roots
      if i != j
        l[j] *= (x - roots[i])/(roots[j] - roots[i])
      end
    end
  end
  return l
end



#=
Collects all solution variables form the DG discretization and averages overlapping nodes
to obtain nodal values for the triangulation solution values.
Edge center node values are the average of the two finder neighbor node values and
the interpolated value on the coarser element using only the edge nodes as interpolation basis.

Also converts the conservative variables into other forms (e.g. primitive variables) if required,
using Trixi functionality (e.g. cons2prim).
=#
function collectSolutionVariables(integrator, active_cells_to_vertices, active_vertices_to_cells, solution_variables)
  n_cells = size(active_cells_to_vertices)[1]
  n_vertices = size(active_vertices_to_cells)[1]

  u_ode = integrator.u
  semi = integrator.p
  #_, equations, solver, _ = Trixi.mesh_equations_solver_cache(semi)
  mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)

  # re-formulate the solution vector into another format
  u = Trixi.wrap_array_native(u_ode, mesh, equations, solver, cache)
  # the solution of a field can then be accessed via data = vec(u[i_var, .., :])
  # (as done in the Trixi solution output)
  # alternatively the vector could be reshaped into a matrix or just split into n=n_variables parts.
  # and then accessed via u[:, i_var] (= e.g. var=1 -> u_ode[1:end/i_var, i_var])
  # Here, to be safe, the same as in Trixi is done


  n_vars = Trixi.nvariables(equations)
  polydeg = Trixi.polydeg(solver)
  nnodes = Trixi.nnodes(solver)
  nodes_per_cell = nnodes*nnodes # 2D quad


  # node indexing:
  # 13 14 15 16 
  #  9 10 11 12
  #  5  6  7  8
  #  1  2  3  4
  # indexing vectors to exctract corner and edge values from the vector of element variable values
  corner_node_positions = [1, nnodes, nodes_per_cell-nnodes+1, nodes_per_cell]
  edge_node_positions = [[1:nnodes:nodes_per_cell...], # left
                         [nnodes:nnodes:nodes_per_cell...], # right
                         [1:nnodes...], # bottom
                         [nodes_per_cell-nnodes+1:nodes_per_cell...]] # top



  # create the coefficients to Lagrange-interpolate the midpoint on the edges
  # (i.e. lagrange basis functions evaluated at x=0 (the edge midpoint in referece space))
  roots = Trixi.LobattoLegendreBasis(polydeg).nodes
  coeffs = lagrangeBasisEvaluation(nnodes, roots, 0.0)


  variables = zeros(n_vertices, n_vars)
  vertex_variables = zeros(8)
  for i_var=1:n_vars
    #data = u_ode[:, i_var] # only works with one variable. Multiple require to split up the vector.
    data = vec(u[i_var, .., :])

    index = 1 # position of the element variables | alternative: in loop: idx_begin = (i_cell-1)*nodes_per_cell + 1
    for i_cell=1:n_cells
      idx_next = index + nodes_per_cell
      element_nodal_values = data[index:idx_next-1]
      index = idx_next

      vertex_variables[1:4] = element_nodal_values[corner_node_positions] # corner values
      for i_edge=1:4 # midpoint values
        if active_cells_to_vertices[i_cell, 4 + i_edge] > 0
          vertex_variables[4 + i_edge] = dot(coeffs, element_nodal_values[edge_node_positions[i_edge]])
        end
      end

      # broadcast to variables matrix
      existing_vertices = active_cells_to_vertices[i_cell, :] .> 0
      vertices = active_cells_to_vertices[i_cell, existing_vertices]
      variables[vertices, i_var] += vertex_variables[existing_vertices]
    end
  end # end for all solution fields (n_vars)

  # average values on the vertices belonging to multiple cells
  variables ./= active_vertices_to_cells[:, end]

  # convert each nodes variables if needed
  if solution_variables != cons2cons
    n_out_vars = size(solution_variables(zeros(n_vars), equations))[1]
    # in place for same output size as input size
    if n_vars == n_out_vars
      for i=1:n_vertices; variables[i,:] = solution_variables(variables[i,:], equations); end

    else
      aux = zeros(n_vertices, n_out_vars)
      for i=1:n_vertices; aux[i,:] = solution_variables(variables[i,:], equations); end
      variables = aux
      # works, but is extremely slow (due to double copy transpose):
      #variables = collect(Matrix(hcat(solution_variables.(eachrow(variables), equations)...))')
    end
  end

  return variables
end








#=

####### assumptions on the mesh | requirements #######

The mesh is a 2D TreeMesh

For all pairs of diagonal neighbors:
One of their shared neighbors active (leaf) cell index is between or lower than the two indices of the diagonal neighbors.
This ensures that the index of a corner vertex is propagated though the cardinal neighbors to the diagonal ones.
For a refined mesh, the same should be true for coarse-fine diagonal neighbors.

There are no interfaces that interface a finer and a coarser cell at once. (Alway max one level difference between cells)
This is auto-satisfied for a tree mesh.

Satisfied by the Trixi framework, TreeMesh and AMR procedure (state Q1 2023):
Corner vertices are nodes of the nodal DG discretization (Gauss-Lobatto quadrature).
A neighboring cell's resolution is at most one level different that the current cell's resolution.
Uniform discretization order (p-conforming).
A reasonable mesh (e.g. no one-cell meshes -> all active leaf cells have a parent)
No inverted coordinate system: x_min=x_left < x_right=x_max (and same for y)
Active cell id (Trixi.leaf_cells(tree)) is is monotonically ascending
An active cell can only have two coarser active neighbors, both connecting to the same corner (true for a TreeMesh)




####### performance #######
Trixi.leaf_cells(tree) used in mapping creation and coordinate creation. -> could get one and pass
use views somewhere?
coordinates: get Trixi.length_at_cell(tree, cell_id) for each level in tree once and store in helper vector
periodic check is done for all cells, irrespective of if the mesh is even periodic
pre-calculate connectivity length when constructing the cell-> vertex mapping or just overestimate size
the triangulation lookup table can be rearranged so that the initial lookup 4D matrix matches 1:16 and can be replaced with a simple function
for repeated execution, a struct storing triangulation tables and the mapping matrices (overestimated size because AMR) would be beneficial
including a dummy row in the mapping matrices at the beginning would eliminate some if > 0 -> just write to i+1 row
major simplifications and performance uplift may be obtained by accounting for to order of active cells in the 2d tree mesh (only need to read left and down and write right and up).
Is the piecewise writing of data to the file slow or fast enough?


=#

