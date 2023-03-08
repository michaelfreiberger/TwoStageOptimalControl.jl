push!(LOAD_PATH,"../src/")
using TwoStageOptimalControl

using Documenter

makedocs(
        sitename="TwoStageOptimalControl.jl",
        modules  = [TwoStageOptimalControl],
        pages=[
                "Home" => "index.md"
               ]
        )

deploydocs(;
        repo = "github.com/michaelfreiberger/TwoStageOptimalControl.jl")

