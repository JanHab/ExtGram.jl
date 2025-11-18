#=
Data collection for 1D tree meshes.
Extends the functionality of the 2D tree triangulator to 1D (although in 1D no triangulation is performed).
=#



#=
Collect Trixi 1D tree mesh data to matrices corresponding to a continuous FE solution
The mapping active_cells_to_vertices is a lambda, to save space.
=#
function collect1dTreeArrays(integrator, solution_variables)
  semi = integrator.p
  u_ode = integrator.u
  mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
  tree = mesh.tree

  n_active_cells = Trixi.count_leaf_cells(tree)
  active_cell_ids = Trixi.leaf_cells(tree)

  basis_nodes = solver.basis.nodes
  n_vertices_per_cell = length(basis_nodes) # Trixi.nnodes(solver)
  # total vertices amount (overlapping vertices combined) = #cell_inner_vertices + #cell_outer_vertices(=n_cells+1)
  n_vertices = (n_vertices_per_cell-1)*n_active_cells+1

  n_vars = Trixi.nvariables(equations)
  u = Trixi.wrap_array_native(u_ode, mesh, equations, solver, cache)

  # i_active_cell -> vertex ids
  # in 1D a simple relation suffices: active_cells_to_vertices(i_cell) = i_cell*(n_vertices_per_cell-1)-3 .+ collect(1:n_vertices_per_cell)
  # use a matrix anyway for ease of use and compatibility (no need to pass along n_active_cells for iteration over cells)
  active_cells_to_vertices = zeros(Int, n_active_cells, n_vertices_per_cell)

  # auxiliary arrays for the solution output
  coordinates = zeros(n_vertices)
  variables = zeros(n_vertices, n_vars)


  # go through all active cells from left to right and compute the nodal coordinates
  for i_cell=1:n_active_cells
    i_cell_global = active_cell_ids[i_cell]
    vertices = (i_cell-1)*(n_vertices_per_cell-1) .+ collect(1:n_vertices_per_cell)
    active_cells_to_vertices[i_cell, :] = vertices
    coordinates[vertices] = basis_nodes*0.5*Trixi.length_at_cell(tree, i_cell_global) .+ Trixi.cell_coordinates(tree, i_cell_global)[1]
  end

  # collect all nodal variable values
  for i_var=1:n_vars
    data = vec(u[i_var, .., :])
    for i_cell=1:n_active_cells
      index = 1 + (i_cell-1)*n_vertices_per_cell # DG -> use a non-overlapping vertex enumeration | alternatively increment indices
      nodal_values = data[index:index+n_vertices_per_cell-1]
      nodal_values[[1,end]] *= 0.5 # average the cell boundary values
      variables[active_cells_to_vertices[i_cell, :], i_var] += nodal_values
    end
  end
  # note: for a periodic boundary, the leftmost and rightmost vertices are NOT averaged
  # undo the 0.5 factor of the leftmost and rightmost vertices, as they are not averaged
  variables[[1,end], :] *= 2

  # convert each nodes variables if needed
  if solution_variables != cons2cons
    n_out_vars = size(solution_variables(zeros(n_vars), equations))[1]
    if n_vars == n_out_vars # in place for same output size as input size
      for i=1:n_vertices; variables[i,:] = solution_variables(variables[i,:], equations); end
    else
      aux = zeros(n_vertices, n_out_vars)
      for i=1:n_vertices; aux[i,:] = solution_variables(variables[i,:], equations); end
      variables = aux
    end
  end

  # create the connectivity for linear FE elements, given by [i, i+1] in each column, since vertices/coordinates are sorted
  # note: just output linear elements in every case for maximum compatibility/support while the outputs remain consistent
  connectivity = collect([1:n_vertices-1 2:n_vertices]') # n_elements = n_vertices-1

  return connectivity, coordinates, variables
end



function collect1DTreeArrays_local(semi, u_ode, solution_variables)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    tree = mesh.tree

    n_active_cells = Trixi.count_leaf_cells(tree)
    active_cell_ids = Trixi.leaf_cells(tree)

    basis_nodes = solver.basis.nodes
    n_vertices_per_cell = length(basis_nodes) # Trixi.nnodes(solver)
    # total vertices amount (overlapping vertices combined) = #cell_inner_vertices + #cell_outer_vertices(=n_cells+1)
    n_vertices = (n_vertices_per_cell-1)*n_active_cells+1

    n_vars = Trixi.nvariables(equations)
    u = Trixi.wrap_array_native(u_ode, mesh, equations, solver, cache)

    # i_active_cell -> vertex ids
    # in 1D a simple relation suffices: active_cells_to_vertices(i_cell) = i_cell*(n_vertices_per_cell-1)-3 .+ collect(1:n_vertices_per_cell)
    # use a matrix anyway for ease of use and compatibility (no need to pass along n_active_cells for iteration over cells)
    active_cells_to_vertices = zeros(Int, n_active_cells, n_vertices_per_cell)

    # auxiliary arrays for the solution output
    coordinates = zeros(n_vertices)
    variables = zeros(n_vertices, n_vars)


    # go through all active cells from left to right and compute the nodal coordinates
    for i_cell=1:n_active_cells
        i_cell_global = active_cell_ids[i_cell]
        vertices = (i_cell-1)*(n_vertices_per_cell-1) .+ collect(1:n_vertices_per_cell)
        active_cells_to_vertices[i_cell, :] = vertices
        coordinates[vertices] = basis_nodes*0.5*Trixi.length_at_cell(tree, i_cell_global) .+ Trixi.cell_coordinates(tree, i_cell_global)[1]
    end

    # collect all nodal variable values
    for i_var=1:n_vars
        data = vec(u[i_var, .., :])
        for i_cell=1:n_active_cells
            index = 1 + (i_cell-1)*n_vertices_per_cell # DG -> use a non-overlapping vertex enumeration | alternatively increment indices
            nodal_values = data[index:index+n_vertices_per_cell-1]
            nodal_values[[1,end]] *= 0.5 # average the cell boundary values
            variables[active_cells_to_vertices[i_cell, :], i_var] += nodal_values
        end
    end
    # note: for a periodic boundary, the leftmost and rightmost vertices are NOT averaged
    # undo the 0.5 factor of the leftmost and rightmost vertices, as they are not averaged
    variables[[1,end], :] *= 2

    # convert each nodes variables if needed
    if solution_variables != cons2cons
        n_out_vars = size(solution_variables(zeros(n_vars), equations))[1]
        if n_vars == n_out_vars # in place for same output size as input size
            for i=1:n_vertices; variables[i,:] = solution_variables(variables[i,:], equations); end
        else
            aux = zeros(n_vertices, n_out_vars)
            for i=1:n_vertices; aux[i,:] = solution_variables(variables[i,:], equations); end
            variables = aux
        end
    end

    # create the connectivity for linear FE elements, given by [i, i+1] in each column, since vertices/coordinates are sorted
    # note: just output linear elements in every case for maximum compatibility/support while the outputs remain consistent
    connectivity = collect([1:n_vertices-1 2:n_vertices]') # n_elements = n_vertices-1

    return connectivity, coordinates, variables
end