module TwoStageOptimalControl

    using ColorSchemes
    using Coverage
    using Dierckx
    using ForwardDiff
    using LaTeXStrings
    using LinearAlgebra
    using Plots
    using Printf
    using Serialization
    using FileIO

    include("AuxiliaryFunction.jl")
    include("LineSearch.jl")
    include("MainFunction.jl")
    include("ModelFunctions.jl")
    include("ParametersVariablesSettings.jl")
    include("ResultsHandling.jl")
    include("StateSolvers.jl")


    export LineSearch
    export TwoStageOptimisation
    export PlotResults
    export SaveResults
    export LoadResults
end