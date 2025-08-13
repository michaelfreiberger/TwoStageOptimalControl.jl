# TwoStageOptimalControl

For general information on the package, the model class it is able to solve, and how to install it, please revisit the [ReadMe page](https://github.com/michaelfreiberger/TwoStageOptimalControl.jl/blob/master/README.md).


## Function Documentation

The following subsections are structured identical to the position of function definitions in the source code and should also act as an orientation for programmers trying to adapt to algorithm for their specific needs.

The most important functions for simple usage are contained in the following three subsections:

- [Main function for usage:](Functions/MainFunction.md) This section contains the function ```TwoStageOptimisation```, which is used to solve a given problem. For examples on how to use it, see [example 1](Examples/Test1.md) and [example 2](Examples/Test2.md).

- [Parameters, Variables, Settings:](Functions/ParametersVariablesSettings.md) In this section all parameters and settings for the solution algorithm are described.

- [Results Handling:](Functions/ResultsHandling.md) In this section all functions are summarised which concern saving, loading and plotting of results.

The remaining subsections contains the documentation for functions used within the solution algorithm:

 - [Auxiliary functions:](Functions/AuxiliaryFunction.md) This section contains functions which are not algorithm specific like integration and smoothing functions.

 - [LineSearch:](Functions/LineSearch.md) This section contains the functions describing the linesearch along the gradient leading to the new guess for the control profiles.

 - [Model Functions:](Functions/ModelFunctions.md) In this section the functions provided by the user are used to define all internal model functions like the Hamiltonian and all equations like the costate dynamics and gradients, which are derived from it.

 - [State solvers:](Functions/StateSolvers.md) In this section all functions for the solution of the state and costate differential equation are defined.

## Examples

As an example often says more then 1000 explanations, please find two examples for the usage of the toolbox below.

 - [Example 1 - Capital accumulation problem](Examples/Test1.md)

 - [Example 2 - Capital accumulation problem with a myopic decision maker](Examples/Test1Myopic.md)

 - [Example 2 - A second capital accumulation problem](Examples/Test2.md)