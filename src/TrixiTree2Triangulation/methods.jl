#=
Entry point of the triangulation callback into the underlying method calls.
Unifies the 2D and 1D conversion (tree flattening and continuous solution creation) and writing functionality.
=#

# 2D
include("converter.jl")
include("writer.jl") # file writing outsourced to separate file
# 1D
include("converter_1d.jl")
include("writer_1d.jl")


#= Main converter function | uses multiple dispatch / specializations to call the corresponding 2D or 1D function =#
function convertMeshAndSolution(integrator, write_options, triangulation_callback)
  mesh, _, _, _ = Trixi.mesh_equations_solver_cache(integrator.p)
  return convertMeshAndSolution(mesh, integrator, write_options, triangulation_callback)
end


#= Main method of this module | uses all above functionality to convert and write the mesh and solution =#
function convertMeshAndSolution(mesh::TreeMesh{2, Trixi.SerialTree{2, T}}, integrator, write_options, triangulation_callback) where {T}
  tree = mesh.tree

  active_cells_to_vertices, n_vertices = generateActiveCellsToVerticesMap(tree)
  active_vertices_to_cells             = generateActiveVerticesToCellsMap(active_cells_to_vertices, n_vertices)
  coordinates                          = collectCoordinates(tree, active_vertices_to_cells)
  connectivity                         = triangulateMesh(active_cells_to_vertices)
  variables                            = collectSolutionVariables(integrator, active_cells_to_vertices,
                                                                  active_vertices_to_cells, triangulation_callback.solution_variables)
  writeFile(write_options, coordinates, variables, connectivity, integrator, triangulation_callback)

  # return some values for verbose callback output (n_vertices, n_triangles)
  return size(coordinates)[1], size(connectivity)[2], integrator.t
end


#= Main method of this module | uses all above functionality to convert and write the 1D mesh and solution =#
function convertMeshAndSolution(mesh::TreeMesh{1, Trixi.SerialTree{1, T}}, integrator, write_options, triangulation_callback) where {T}
  connectivity, coordinates, variables = collect1dTreeArrays(integrator, triangulation_callback.solution_variables)
  writeFile(write_options, coordinates, variables, connectivity, integrator, triangulation_callback)
  return size(coordinates)[1], size(connectivity)[2], integrator.t
end
