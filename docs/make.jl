import Pkg
using Logging
Pkg.activate(@__DIR__)
Pkg.develop(; path=joinpath(@__DIR__, ".."))
Pkg.instantiate()

ENV["GKSwstype"] = "100"

using Documenter
using ExtGram

DocMeta.setdocmeta!(ExtGram, :DocTestSetup, :(using ExtGram); recursive=true)

makedocs(
    modules = [ExtGram],
    sitename = "ExtGram.jl",
    remotes = nothing,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        repolink = "https://github.com/JanHab/ExtGram.jl",  
        edit_link = "main"
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "1D1D - Shock Tube Example" => "examples/1D1D_ShockTube.md",
            "1D1D - Shock Structure Example" => "examples/1D1D_ShockStructure.md",
            "1D3V - Shock Tube Example" => "examples/1D3V_ShockTube.md",
        ],
        "API Reference" => "reference.md",
    ],
    clean = true,
    # strict = true,
)
