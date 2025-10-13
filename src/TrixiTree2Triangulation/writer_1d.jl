#=
File writer specializations for the 1D case.
Specialization is performed via multiple dispatch on coordinates::Vector (Matrix in 2D case).
(except for the hdf5 writer, which does not need changes)
In all cases, the changes are minimal (fancier code-reuse architecture is avoided to save time).
=#

#= tab separated plain text data file of 1D solution for use with e.g. gnuplot =#
function write_tsv(coordinates::Vector, variables, connectivity, write_options, metadata, integrator, callback)
  _, _, append, first_write, _, = write_options
  filepath, timestep, time, equations_name, variable_names, n_vertices, _, n_variables = metadata
  file_config = "w"; if (append && !first_write); file_config = "a"; end


  open(filepath, file_config) do file
    # header | only write the file header once in append mode
    if !append || first_write
      write(file, "# \"Trixi 1D tree mesh continuous solution (Equations: $equations_name)\"\n")
      write(file, "# info=\"" * callback.info_string * "\"\n")
      write(file, "# \"X\"")
      for i_var=1:n_variables; write(file, ", \"", variable_names[i_var], "\""); end; write(file, "\n");
    end

    # zone header
    write(file, "# timestep=$timestep, points=$n_vertices, t=$time\n")

    # coordinates and variables
    for i_vertex=1:n_vertices
      write(file, @sprintf "%.14E" coordinates[i_vertex])
      for i_var=1:n_variables
        write(file, " ", @sprintf "%.14E" variables[i_vertex, i_var])
      end
      write(file, "\n")
    end
    # leave a blank line at the end (doubles as block separator for gnuplot)
    write(file, "\n")
  end # for open file do
end





#= 1D version of the tecplot .dat format writer =#
function write_dat(coordinates::Vector, variables, connectivity, write_options, metadata, integrator)
  _, _, append, first_write, _, = write_options
  filepath, timestep, time, equations_name, variable_names, n_vertices, n_cells, n_variables = metadata
  file_config = "w"; if (append && !first_write); file_config = "a"; end


  open(filepath, file_config) do file
    # header | only write the file header once in append mode
    if !append || first_write
      write(file, "TITLE=\"Trixi 1D tree mesh continuous solution in Tecplot (Equations: $equations_name)\"\n")
      write(file, "VARIABLES = X")
      for i_var=1:n_variables; write(file, ", \"", variable_names[i_var], "\""); end; write(file, "\n");
    end

    # zone header
    write(file, "ZONE T=\"ts$timestep\", NODES=$n_vertices, ELEMENTS=$n_cells, DATAPACKING=POINT, ZONETYPE=FELINESEG, SOLUTIONTIME=$time\n")

    # coordinates and variables
    for i_vertex=1:n_vertices
      write(file, @sprintf "%.14E" coordinates[i_vertex])
      for i_var=1:n_variables
        write(file, " ", @sprintf "%.14E" variables[i_vertex, i_var])
      end
      write(file, "\n") # point data-packing: end line after all info on a point
    end

    # connectivity
    write(file, "\n")
    for i_elem=1:n_cells
      # vertices/coordinates and elements are sorted
      write(file, "$(i_elem) $(i_elem+1)\n")
      #write(file, "$(connectivity[1, i_elem]) $(connectivity[2, i_elem])\n")
    end
  end # for open file do

end



#= 1D version of the tecplot .szplt format writer =#
function write_szplt(coordinates::Vector, variables, connectivity, write_options, metadata, integrator)
  _, _, append, first_write, file_handle = write_options
  filepath, timestep, time, equations_name, variable_names, n_vertices, n_cells, n_variables = metadata

  #= open file and create header | only done once in append mode =#
  if !append || first_write
    var_string = "X";
    for i_var=1:n_variables; var_string = string(var_string, " ", variable_names[i_var]); end;
    data_set_title = "Trixi 1D tree mesh continuous solution in Tecplot (Equations: $equations_name)"

    # file format .szplt (=1), file type full (=0), default data type double (=2)
    ccall((:tecFileWriterOpen, libtecio), Cint, (Cstring, Cstring, Cstring, Cint, Cint, Cint, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}),
                                                 filepath, data_set_title, var_string, 1, 0, 2, C_NULL, file_handle)
  end


  #= Create a new zone for the solution data and mesh of the current time step =#
  value_locations = ones(Cint, n_variables+2); # same for each variable, =1 means nodal values
  value_locations_ptr = pointer(value_locations)
  zone_handle = Ref{Cint}() # index and handle of the zone | note: access a reference via zone_handle[]
  zone_title = "ts$timestep"
  # zone type FE line segment (=1)
  ccall((:tecZoneCreateFE, libtecio), Cint, (Ptr{Cvoid}, Cstring, Cint, Int64, Int64,
                                             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                                             Cint, Int64, Cint, Ref{Cint}),
                                             file_handle[], #file_handle[] dereference the file_handle (expects *void)
                                             zone_title, 1, n_vertices, n_cells,
                                             C_NULL, C_NULL, value_locations_ptr, C_NULL,
                                             0, 0, 0, zone_handle)
  # set solution time
  ccall((:tecZoneSetUnsteadyOptions, libtecio), Cint, (Ptr{Cvoid}, Cint, Float64, Cint), file_handle[], zone_handle[], time, 1)


  #= Write the coordinates, solution data and connectivity =#
  ccall((:tecZoneVarWriteDoubleValues, libtecio), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Int64, Ptr{Float64}),
                                                        file_handle[], zone_handle[], 1, 1, n_vertices, pointer(coordinates));
  for i=1:n_variables
    j = i+1 # offset variable counter by X-coord
    ccall((:tecZoneVarWriteDoubleValues, libtecio), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Int64, Ptr{Float64}),
                                                          file_handle[], zone_handle[], j, 1, n_vertices, pointer(variables[:,i]));
  end
  # column major storage yields the correct column-wise (one element after the other) access to the connectivity matrix
  ccall((:tecZoneNodeMapWrite64, libtecio), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Int64, Ptr{Int64}),
                                                  file_handle[], zone_handle[], 1, 1, 2*n_cells, pointer(connectivity));

  # flush the file, storing all currently collected data into the file
  ccall((:tecFileWriterFlush, libtecio), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}), file_handle[], 1, zone_handle)


  #= close the file if we do not append next step or if this is the final step =#
  if !append || Trixi.isfinished(integrator)
    ccall((:tecFileWriterClose, libtecio), Cint, (Ptr{Cvoid},), file_handle)
  end
end



#= 1D version of the vtu writer =#
function write_vtu(coordinates::Vector, variables, connectivity, write_options, metadata, integrator)
  filepath, timestep, time, equations_name, variable_names, _, _, n_variables = metadata

  # the use of VTL_LINE here is the only difference to the 2D version
  cells = map((c)-> MeshCell(VTKCellTypes.VTK_LINE, c), eachcol(connectivity));

  vtk_grid(filepath, coordinates', cells) do file
    for i=1:n_variables; file[variable_names[i]] = variables[:, i]; end
    file["time"] = time
    file["equations"] = string(equations_name)
  end
end
