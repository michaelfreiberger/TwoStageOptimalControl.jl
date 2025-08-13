# Capital accumulation problem for a myopic planer

We can also use the solver to solve the capital accumulation problem from [Example 1](Test1.md) for a myopic planer. Myopic planer means that the planer is not aware of the risk of a shock occuring and behaves accordingly. However, once the shock has occured the planer adjusts its behaviour according to the new system.

## First stage solution

We start by solving the first stage problem without shock. This can be easily achieved by setting the dimension of the state and control variables of the second stage equal to zero.

```julia-repl
julia> MyPara["nCon_dist"] = 0
julia> MyPara["nStat_dist"] = 0
```

Note, that we could theoretically also reduce the dimensions of the state variables by 1, as we do not need to model the survival probability in this scenario. However, calculating and following the dynamics of one additional state variable does not significantly increase the computation effort and first and foremost has no effect on the optimal solution.

After defining the functions for the first stage as in [Example 1](Test1.md)
```julia-repl
julia> U(Con, Stat, t::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
julia> f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1]]
julia> S1(Stat, Para::Dict) = 0
```
we can again define some intial guesses for the control profiles. This is optional, but if we want to provide a starting guess for the first stage solution, we also have to provide an intial guess for the second stage variables. However, note that the last dimension of ```Results["Con_dist"]$``` is equal to zero.

```julia-repl
julia> Results = Dict()
julia> Results["Con"] = 4.0*ones(1,100,1)
julia> Results["Con_dist"] = 4.0*ones(100,100,0)
```

Finally, we can call the main function with just providing the first stage functions.
```julia-repl
julia> Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand1 = U, 
                                StateDynamic1 = f1,
                                SalvageFunction1=S1)
```


## Second stage solution

For the optimal reactions after the shock, we can follow the strategy below. First, we again define the system dynamics of the second stage as in [Example 1](Test1.md).

```julia-repl
julia> Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
julia> f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1], 0]
julia> g(Con,Stat,t::Float64, Para::Dict) = [(1-Para["gamma"])*Stat[1],Para["eta"]*Stat[2]]
```
Furthermore, we need to adjust the parameters of the dimensions of the second stage to their actual values:
```julia-repl
julia> MyPara["nCon_dist"] = 1
julia> MyPara["nStat_dist"] = 2
```

We can then use the option to set the parameters "FixedControl1" equal to a list of all dimensions of the model (in this case just 1)

```julia-repl
julia> Para2["FixedControl1"] = [1]
```

If we now provide the optimal solution of the first-stage-only solution as an initial guess to the solver, this solution will not be changed in the process. We just need to define the initial guess for the second stage solution:

```julia-repl
julia> Results["Con_dist"] = 4.0*ones(size(Results["Con"])[2],size(Results["Con"])[2],1)
```

If we consequently call the main optimisation function, the second stage gets solved for all possible shock timings simultaneously.
```julia-repl
julia> Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand1 = U, 
                                AggregationIntegrand2 = Q,
                                StateDynamic1 = f1,
                                StateDynamic2 = f2, 
                                Shock12 = g,
                                SalvageFunction1=S1,
                                SalvageFunction2=S2)
```

And again, we can save and plot the results.

```julia-repl
julia> SaveResults(Results,"CapitalAccumulationMyopic")
julia> PlotResults(Results,SavePlot = true)
```
