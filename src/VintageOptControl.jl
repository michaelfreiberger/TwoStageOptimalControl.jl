module TwoStageOptControl

    using LinearAlgebra
    using ForwardDiff
    using Interpolations

    using Plots
    using ColorSchemes
    using LaTeXStrings
    using Printf

    using Serialization
    using Dierckx
    using BenchmarkTools
    using Profile
    using Coverage

    include("AuxiliaryFunction.jl")
    include("LineSearch.jl")
    include("MainFunction.jl")
    include("ModelFunctions.jl")
    include("ParametersVariablesSettings.jl")
    include("ResultsHandling.jl")
    include("StateSolvers.jl")


    export LineSearch
    export VintageOptimisation
    export PlotResults

end
