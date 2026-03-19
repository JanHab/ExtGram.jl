#=
Wraps the triangulation and output functionality implemented in converter.jl in a callback object.

Modeled after and mirrors sections from SaveSolutionCallback{SolutionVariables} in Trixi: save_solution.jl

Useful Trixi files:
https://github.com/trixi-framework/Trixi.jl/blob/main/src/callbacks_step/save_solution.jl
https://github.com/trixi-framework/Trixi.jl/blob/main/src/callbacks_step/save_solution_dg.jl

https://github.com/trixi-framework/Trixi.jl/blob/24a03e360ce730aa0fc8b6ddbe509f1a64ef7351/examples/tree_2d_dgsem/elixir_advection_callbacks.jl
=#


"""
Callback struct that triangulates and writes the 2D tree mesh and the solution variables.
Can also be used for 1D Tree meshes, in which case a continuous piecewise linear solution is constructed.

#Configuration:

-`solution_variables`: The set of solution variables to write (via Trixi solution variables function, such as 'cons2prim')
-`interval`: Save every interval time-steps. An interval <= 0 means do not save interval-based.

-`time_interval`: Save at every simulation time multiple of the time interval.
A time interval <= 0 disables this functionality. Cannot be used together with interval based saving.
The callback guarantees that the interval multiples are hit by reducing the time-step size if needed.

-`save_initial_solution`, `save_final_solution`, `output_directory`

-`file_format`: supported file-formats are "dat", "szplt", "vtu" and "h5".
  '.dat' is the Tecplot Tabular data file written in ascii, readable by both Tecplot and Paraview.
  '.szplt' is the Tecplot proprietary format only readable by Tecplot.
    Use of szplt writing requires the LD_LIBRARY_PATH to contain the libtecio.so shared library file.
    Alternatively libtecio.so can also be present in the working directory.
  '.vtu' is the native Paraview format for an unstructured grid.
  '.h5' is a dump of the relevant data matrices into a hdf5 file, for more direct reading with e.g. Mathematica.
  '.tsv' 1D only file format containing the vertex data
  'trixi' use Trixi.SaveSolutionCallback to generate the output, instead of this triangulation callback.

-`append_solution`: If true, write all solution data into one file.
  Only the dat, szplt, h5 and tsv file formats support the append option. vtu will simply ignore the append flag.
  Note that Tecplot only recognizes time evolutions when they are in one file, while Paraview is the exact opposite.

-`verbose`: print a short message to console every time output is written

-`clear_out_dir`: delete the contents of the output directory at initialization
"""
mutable struct SaveTriangulationCallback{SolutionVariables}
  # configuration
  solution_variables::SolutionVariables
  interval::Int
  time_interval::Float64
  save_initial_solution::Bool
  save_final_solution::Bool
  output_directory::String
  file_format::String # string simpler than enum or class
  append_solution::Bool
  verbose::Bool
  clear_out_dir::Bool
  file_basename::String
  info_string::String # tsv file format only
  # helpers
  first_write::Bool # set to true in initializer
  file_handle::Ref{Ptr{Cvoid}} # **void handle for a tecio file
end


# this method is called when the callback is activated
function (triangulation_callback::SaveTriangulationCallback)(integrator)

  Trixi.@trixi_timeit Trixi.timer() "triangulation & output" begin
    write_options = (triangulation_callback.output_directory,
                     triangulation_callback.file_format,
                     triangulation_callback.append_solution,
                     triangulation_callback.first_write,
                     triangulation_callback.file_handle)
    triangulation_callback.first_write = false
    nv, nt, time = convertMeshAndSolution(integrator, write_options, triangulation_callback)

    if triangulation_callback.verbose
      println("Triangulation callback: written output (time = ", round(time; digits = 4),", #vertices = ", nv, ", #triangles = ", nt, ")")
    end
  end

  # avoid re-evaluation possible FSAL stages (First Same As Last, Dormand–Prince method)
  u_modified!(integrator, false)
  return nothing
end



# determine whether the discrete callback should be called
function (triangulation_callback::SaveTriangulationCallback)(u, t, integrator)
  @unpack interval, time_interval, save_final_solution = triangulation_callback
  #interval = triangulation_callback.interval
  #save_final_solution = triangulation_callback.save_final_solution
  
  # only count the accepted time-steps (rejections due to error-based time-step control)
  # thus some fancy conditions are used (see Trixi: save_solution.jl)
  interval_hit = (interval > 0 # interval saving enabled
                  && (integrator.stats.naccept % interval == 0) # interval hit
                  && (integrator.stats.naccept > 0)) # skip initial step(s)
  is_finished = (save_final_solution && Trixi.isfinished(integrator))

  return interval_hit || is_finished
end


# helper to create or clear the output directory (if needed) during initialization
function manage_output_directory(output_directory, clear_out_dir, verbose)
  if Trixi.mpi_isroot()
    if clear_out_dir
      if verbose; println("Clearing output directory \"", output_directory, "\""); end
      rm(output_directory, recursive=true, force=false)
    end
    mkpath(output_directory)
  end
end


# initialization of the triangulation callback and activation if save_initial_solution
function Trixi.initialize_save_cb!(triangulation_callback::SaveTriangulationCallback, u, t, integrator)

  # check wether the specified file format is supported
  supported = ["dat", "szplt", "h5", "vtu", "tsv", "trixi"]
  if !(triangulation_callback.file_format in supported)
    @error string("Triangulation Callback: specified file format ",
                  triangulation_callback.file_format, " is not supported. Allowed values are ",
                  string("'", join(supported, "' '"), "'"))
  end
  # check tsv only for 1D
  if triangulation_callback.file_format == "tsv" && Trixi.ndims(integrator.p) != 1
    @error string("Triangulation Callback: tsv file format only supported for 1D tree meshes")
  end

  # for vtu format, disable appending (unsupported, but default = true)
  if triangulation_callback.file_format == "vtu"; triangulation_callback.append_solution = false; end

  if triangulation_callback.info_string != "" && triangulation_callback.file_format!="tsv"
    @warn "Triangulation callback: 'info' keyword argument only supported for tsv files."
  end

  # create the directory if needed and write initial solution if configured
  manage_output_directory(triangulation_callback.output_directory,
                          triangulation_callback.clear_out_dir,
                          triangulation_callback.verbose)
  if triangulation_callback.save_initial_solution
    triangulation_callback(integrator)
  end

  return nothing
end


# Create a 'DiscreteCallback' for the callback set
function SaveTriangulationCallback(; solution_variables=Trixi.cons2cons,
                                     interval=0,
                                     time_interval=0.0,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     output_directory="out",
                                     file_format="dat",
                                     append_solution=true,
                                     verbose=true,
                                     clear_out_dir=false,
                                     name="solution",
                                     info="")

  # if configured, use the standard Trixi output callback instead of the triangulation
  if file_format=="trixi"
    # clear the output directory here, since SaveSolutionCallback cannot
    manage_output_directory(output_directory, clear_out_dir, verbose)
    # switch to the purely the Trixi callback
    if verbose
      println("Triangulation callback: switching to Trixi.SaveSolutionCallback. No future verbose output.")
    end
    return triangulation_callback = Trixi.SaveSolutionCallback(interval=interval,
                                                        dt=time_interval,
                                                        save_initial_solution=save_initial_solution,
                                                        save_final_solution=save_final_solution,
                                                        output_directory=output_directory,
                                                        solution_variables=solution_variables)
  end

  triangulation_callback = SaveTriangulationCallback(solution_variables,
                                                     interval,
                                                     time_interval,
                                                     save_initial_solution,
                                                     save_final_solution,
                                                     output_directory,
                                                     file_format,
                                                     append_solution,
                                                     verbose,
                                                     clear_out_dir,
                                                     name,
                                                     info,
                                                     true, # first_write=true
                                                     Ref{Ptr{Cvoid}}()) # null ptr for tecio file handle
  #
  if time_interval > 0.0
    if interval > 0
      @error "Triangulation Callback: interval and time interval based saving at once not supported, use two separate callbacks"
    end
    return PeriodicCallback(triangulation_callback, time_interval,
                            save_positions = (false, false),
                            initialize = Trixi.initialize_save_cb!, # recursive multiple dispatch calls the implementation above
                            final_affect = save_final_solution)
  else
    # first triangulation_callback is the struct-method() overload for the condition, the second one the actual callback
    return DiscreteCallback(triangulation_callback, triangulation_callback,
                            save_positions=(false,false),
                            initialize = Trixi.initialize_save_cb!)
  end
end
