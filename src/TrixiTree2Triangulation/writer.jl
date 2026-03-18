#=
File output of the triangulation matrices created in converter.jl
Supports functions for different file formats.
Used in converter.jl convertMeshAndSolution()
=#

#=
Write the triangulated mesh and solution into a file.
Switch function for the different file format writers.

A more elegant approach would be multiple dispatch specialization of this function
based on file-format object specializations or enum values.
This would also make callback file format configuration smoother.
=#
function writeFile(write_options, coordinates, variables, connectivity, integrator, callback)
  ### Do some prep that is needed for all writers here ###
  directory, format, append, _, _ = write_options
  file_basename = callback.file_basename

  # unpack integrator and solution
  timestep = integrator.stats.naccept
  time = integrator.t
  semi = integrator.p
  _, equations, _, _ = Trixi.mesh_equations_solver_cache(semi)
  equations_name = Base.typename(typeof(equations)).wrapper
  variable_names = Trixi.varnames(callback.solution_variables, equations)
  n_vertices = size(coordinates)[1]
  n_triangles = size(connectivity)[2]
  n_variables = size(variables)[2]

  # create filename | filename for single file solution is just "solution.format" while for non-append, we suffix the time-step
  filename = string(file_basename, ".", format)
  if !append
    format_string = string(file_basename,"_%06d.")
    filename = string(Printf.format(Printf.Format(format_string), timestep), format)
  end
  filepath = joinpath(directory, filename)

  metadata = (filepath, timestep, time, equations_name, variable_names, n_vertices, n_triangles, n_variables)

  if format == "dat"       write_dat(coordinates, variables, connectivity, write_options, metadata, integrator);
  elseif format == "szplt" write_szplt(coordinates, variables, connectivity, write_options, metadata, integrator);
  elseif format == "h5"    write_hdf5(coordinates, variables, connectivity, write_options, metadata, integrator);
  elseif format == "vtu"   write_vtu(coordinates, variables, connectivity, write_options, metadata, integrator);
  elseif format == "tsv"   write_tsv(coordinates, variables, connectivity, write_options, metadata, integrator, callback);
  else; @warn "Unsupported file format, skipping output." # re-warn (callback shows error if unsupported)
  end
end








#=
Write the triangulated mesh and solution into an ascii .dat tecplot data file.
Uses DATAPACKING=POINT.
Currently writes out 14 digits of the coordinates and variables.
=#
function write_dat(coordinates, variables, connectivity, write_options, metadata, integrator)
  # unpack write_options
  directory, _, append, first_write, _, = write_options
  file_config = "w"; if (append && !first_write); file_config = "a"; end

  filepath, timestep, time, equations_name, variable_names, n_vertices, n_triangles, n_variables = metadata

  

  open(filepath, file_config) do file
    # header | only write the file header once in append mode
    if !append || first_write
      write(file, "TITLE=\"Trixi triangulated to Tecplot (Equations: $equations_name)\"\n")
      write(file, "VARIABLES = X, Y")
      for i_var=1:n_variables; write(file, ", \"", variable_names[i_var], "\""); end; write(file, "\n");
    end

    # zone header
    write(file, "ZONE T=\"ts$timestep\", NODES=$n_vertices, ELEMENTS=$n_triangles, DATAPACKING=POINT, ZONETYPE=FETRIANGLE, SOLUTIONTIME=$time\n")

    # coordinates and variables
    for i_vertex=1:n_vertices
      write(file, @sprintf "%.14E" coordinates[i_vertex, 1])
      write(file, " ", @sprintf "%.14E" coordinates[i_vertex, 2])
      for i_var=1:n_variables
        write(file, " ", @sprintf "%.14E" variables[i_vertex, i_var])
      end
      write(file, "\n") # point data-packing: end line after all info on a point
    end

    # connectivity
    write(file, "\n")
    for i_elem=1:n_triangles
      write(file, "$(connectivity[1, i_elem]) $(connectivity[2, i_elem]) $(connectivity[3, i_elem])\n")
    end
  end # for open file do

end






const libtecio = "libtecio.so" # library file name | local definition forbidden
#=
Write the triangulated mesh and solution into a binary .szplt tecplot data file using the TecIO.c library.
Currently writes double precision (Float64 and Int64) data and connectivity information.

Useful references concerning the use of TecIO in Julia:
Tecplot360EX Data Format Guide (Tecplot 360 EX 2022 Release 2):
https://tecplot.azureedge.net/products/360/current/360_data_format_guide.pdf
tecio flushpartitioned example: <tecplot_folder>/utils/tecio/examples/flushpartitioned/flushpartitioned.cpp
https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/
https://stackoverflow.com/questions/40140699/the-proper-way-to-declare-c-void-pointers-in-julia
=#
function write_szplt(coordinates, variables, connectivity, write_options, metadata, integrator)
  # unpack write_options
  directory, _, append, first_write, file_handle = write_options

  filepath, timestep, time, equations_name, variable_names, n_vertices, n_triangles, n_variables = metadata

  #= open file and create header | only done once in append mode =#
  if !append || first_write
    var_string = "X Y";
    for i_var=1:n_variables; var_string = string(var_string, " ", variable_names[i_var]); end;
    data_set_title = "Trixi triangulated to Tecplot (Equations: $equations_name)"

    # file format .szplt (=1), file type full (=0), default data type double (=2)
    ccall((:tecFileWriterOpen, libtecio), Cint, (Cstring, Cstring, Cstring, Cint, Cint, Cint, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}),
                                                 filepath, data_set_title, var_string, 1, 0, 2, C_NULL, file_handle)
    #
    # enable diagnostic output
    #ccall((:tecFileSetDiagnosticsLevel, libtecio), Cint, (Ptr{Cvoid}, Cint), file_handle[], 1)
  end


  #= Create a new zone for the solution data and mesh of the current time step =#
  value_locations = ones(Cint, n_variables+2); # same for each variable, =1 means nodal values
  value_locations_ptr = pointer(value_locations)
  zone_handle = Ref{Cint}() # index and handle of the zone | note: access a reference via zone_handle[]
  zone_title = "ts$timestep"
  # zone type FE triangle (=2)
  ccall((:tecZoneCreateFE, libtecio), Cint, (Ptr{Cvoid}, Cstring, Cint, Int64, Int64,
                                             Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                                             Cint, Int64, Cint, Ref{Cint}),
                                             file_handle[], #file_handle[] dereference the file_handle (expects *void)
                                             zone_title, 2, n_vertices, n_triangles,
                                             C_NULL, C_NULL, value_locations_ptr, C_NULL,
                                             0, 0, 0, zone_handle)
  # set solution time
  ccall((:tecZoneSetUnsteadyOptions, libtecio), Cint, (Ptr{Cvoid}, Cint, Float64, Cint), file_handle[], zone_handle[], time, 1)


  #= Write the coordinates, solution data and connectivity =#
  for i=1:2
    ccall((:tecZoneVarWriteDoubleValues, libtecio), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Int64, Ptr{Float64}),
                                                          file_handle[], zone_handle[], i, 1, n_vertices, pointer(coordinates[:,i]));
  end
  for i=1:n_variables
    j = i+2 # offset variable counter by X,Y
    ccall((:tecZoneVarWriteDoubleValues, libtecio), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Int64, Ptr{Float64}),
                                                          file_handle[], zone_handle[], j, 1, n_vertices, pointer(variables[:,i]));
  end
  # we rely on column major storage and each column representing one triangle for the correct linear access to the matrix
  ccall((:tecZoneNodeMapWrite64, libtecio), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Int64, Ptr{Int64}),
                                                  file_handle[], zone_handle[], 1, 1, 3*n_triangles, pointer(connectivity));


  # flush the file, storing all currently collected data into the file
  ccall((:tecFileWriterFlush, libtecio), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}), file_handle[], 1, zone_handle)


  #= close the file if we do not append next step or if this is the final step =#
  if !append || Trixi.isfinished(integrator)
    ccall((:tecFileWriterClose, libtecio), Cint, (Ptr{Cvoid},), file_handle)
  end
end



#=
Write the triangulated mesh and solution into a hdf5 .h5 file.
The output is simply an organized dump of all relevant arrays/matrices.
For easy import into other software, such as Mathematica.

This function is the same for 2D and 1D data.
=#
function write_hdf5(coordinates, variables, connectivity, write_options, metadata, integrator)
  # unpack write_options
  _, _, append, first_write, _, = write_options
  file_config = "w"; if (append && !first_write); file_config = "cw"; end

  filepath, timestep, time, equations_name, variable_names, n_vertices, n_triangles, n_variables = metadata


  h5open(filepath, file_config) do file
    if first_write; attributes(file)["equations"] = string(equations_name); end

    # create a group for the current time-step
    g = create_group(file, "ts$timestep")
    # write vertices and connectivity as datasets
    g["vertices"] = coordinates
    g["connectivity"] = connectivity
  
    # write the individual variables into datasets with corresponding names
    gg = create_group(g, "variables")
    for i=1:n_variables; gg[variable_names[i]] = variables[:, i]; end
    close(gg)

    # write the metadata (i.e. solution time) into an attribute of the time-step group
    attributes(g)["time"] = time
    close(g)
  end
end




#=
Write the triangulated mesh and solution into a Paraview .vtu file using the WriteVTK julia package.
=#
function write_vtu(coordinates, variables, connectivity, write_options, metadata, integrator)

  filepath, timestep, time, equations_name, variable_names, n_vertices, n_triangles, n_variables = metadata

  # map all connectivity triangles to VTK cells
  cells = map((c)-> MeshCell(VTKCellTypes.VTK_TRIANGLE, c), eachcol(connectivity));

  vtk_grid(filepath, coordinates', cells) do file
    for i=1:n_variables; file[variable_names[i]] = variables[:, i]; end
    file["time"] = time
    file["equations"] = string(equations_name)
  end
end
