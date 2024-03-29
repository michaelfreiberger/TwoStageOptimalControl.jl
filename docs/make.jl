#push!(LOAD_PATH,"../src/")
using Documenter, TwoStageOptimalControl

makedocs(
        sitename="TwoStageOptimalControl.jl",
        modules  = [TwoStageOptimalControl],
        pages=[
        "Home" => "index.md",
        "Main Functions" => "Functions/MainFunction.md",
        "Parameters, Variables, Settings" => "Functions/ParametersVariablesSettings.md",
        "Results analysis" => "Functions/ResultsHandling.md",
        "Auxiliary functions" => "Functions/AuxiliaryFunction.md",
        "LineSearch" => "Functions/LineSearch.md",
        "Model Functions" => "Functions/ModelFunctions.md",
        "State Solvers" => "Functions/StateSolvers.md",
        "Example 1" => "Examples/Test1.md",
        "Example 2" => "Examples/Test3.md",
               ]
        )

deploydocs(;
       repo = "github.com/michaelfreiberger/TwoStageOptimalControl.jl.git",
)

 