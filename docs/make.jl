import Pkg
using Logging
Pkg.activate(@__DIR__)
# try
#     Pkg.Registry.add("General")
# catch err
#     @info "General registry already available" err
# end
# try
#     Pkg.Registry.update()
# catch err
#     @warn "Registry update failed" err
# end
Pkg.develop(; path=joinpath(@__DIR__, ".."))
Pkg.instantiate()

ENV["GKSwstype"] = "100"

using Documenter
using HyQMOM

DocMeta.setdocmeta!(HyQMOM, :DocTestSetup, :(using HyQMOM); recursive=true)

makedocs(
    modules = [HyQMOM],
    sitename = "HyQMOM.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true"),
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Shock Tube Example" => "examples/main.md",
        ],
        "API Reference" => "reference.md",
    ],
    clean = true,
    # strict = true,
)
