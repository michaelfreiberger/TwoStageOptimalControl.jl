# Capital accumulation problem

Consider the following capital accumulation problem, where a stochastic shock destroys are certain share of the production capital. The problem consists of the following variables:

- Production capital $K_1(t)$ before the shock (first stage), $K_2(t,s)$ after the shock at time $s$ (second stage).
- Investment in production capital in the first, $I_1(t)$, and second, $I_2(t,s)$, stage.
- Consumption $c_1(t)$ and $c_2(t,s)$ in the first and second stage. Consumption is thereby the residual of production minus the cost for investment in production capital (which is assumed to be quadratic)

```math
\begin{aligned}
c_1(t) &= K_1(t)^{0.5} - I_1(t)^2 \\
c_2(t,s) &= K_2(t,s)^{0.5} - I_2(t,s)^2
\end{aligned}
```

- The shock has an arrival rate of $\eta$ with the time horizon being scaled down to $T=1.0$. At the point of the shock, we assume a share $\gamma$ of production capital is being destroyed.

The social planer now maximizes the aggregated consumption over the whole time horizon. The whole problem writes as follows:

```math
\begin{aligned}
\max_{I_1(t),I_2(t,s)} &\int_0^T e^{\rho t}\cdot\left\{ Z_1(t)\cdot \Big[ K_1(t)^{0.5} - I_1(t)^2\Big] + Q(t) \right\} dt \\
s.t.\hspace{1cm} \dot{K_1}(t) &= I_1(t) - \delta K_1(t) &&,\qquad K_1(0) = K_1^0 \\
\dot{Z_1}(t) &= - \eta\cdot Z_1(t)                      &&,\qquad Z_1(0) = 1 \\
\dot{K_2}(t,s) &= I_2(t,s) - \delta K_2(t,s)            &&,\qquad K_2(s,s) = (1-\gamma)\cdot K_1(s) \\
\dot{Z_2}(t,s) &= 0                                     &&,\qquad Z_2(s,s) = \eta\cdot Z_1(s) \\
Q(t) &= \int_0^t Z_2(t,s)\cdot \Big[ K_2(t,s)^{0.5} - I_2(t,s)^2\Big] ds
\end{aligned}
```

# Solution in Julia

To solve this problem with our Toolbox, we can use the main function ```TwoStageOptimisation```.

```julia-repl
julia> Results = TwoStageOptimisation(UserParameters = MyPara,
                                ObjectiveIntegrand = U, 
                                AggregationFunction = Q,
                                StateDynamic_1 = f1,
                                StateDynamic_2 = f2, 
                                Shock = g,
                                SalvageFunction_1=S1,
                                SalvageFunction_2=S2)
```

But first we need to define the functions ```U, Q, f1, f2, g, S1``` and ```S2``` as well as the dictionary ```MyPara```.

Let us start with the dictionary.

## MyPara

The dictionary ```MyPara``` contains all relevant parameters of the model. We can initialize an empty dictionary by
```julia-repl
julia> MyPara = Dict()
```
```MyPara``` contains information on

#### Model specific parameters

like $\rho$, $\eta$ and $\gamma$ in this example;
```julia-repl
julia> MyPara["T"] = 1
julia> MyPara["rho"] = 0.05
julia> MyPara["delta"] = 0.1
julia> MyPara["eta"] = 0.5
julia> MyPara["gamma"] = 0.2
```
#### Model dimensions of controls and states

```nCon``` is the number of control variables in the first stage, ```nCon_dist``` is the number of control variables in the second stage. ```nStat``` and ```nStat_dist``` are the dimensions of state variables before and after the shock.

```julia-repl
julia> MyPara["nCon"] = 1
julia> MyPara["nCon_dist"] = 1
julia> MyPara["nStat"] = 2
julia> MyPara["nStat_dist"] = 2
```
#### Control boundaries and initial state values

```julia-repl
julia> MyPara["ConMin"] = [0.0]
julia> MyPara["Con_distMin"] = [0.0]
julia> MyPara["InitStat"] = [0.1,1.0]
```

#### Algorithm specifying parameters

```julia-repl
julia> MyPara["hstep"] = 0.05
julia> MyPara["InitLineStep"] = 1e-5
julia> MyPara["UpperLineStep"] = 1e-1
julia> MyPara["hLowBound"] = 0.01
julia> MyPara["PlotResultsIntermediateFrequency"] = 30
```

Information on the parameters base values, which are used if not specified differently, can be found in the section [Parameters, Variables, Settings](../Functions/ParametersVariablesSettings.md)

## Model functions

Now we need to specify all the functions used in the model. Any function, that is not defined is assumed to be the zero-function.

* ```U(Con, Stat, t::Float64, Para::Dict)``` covers the utility function in the first stage (without the aggregated part $Q(t)$ of the second stage). In our example this is 
```julia-repl
julia> U(Con, Stat, t::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
```

* ```Q(Con,Stat, t::Float64,s::Float64, Para::Dict)``` contains the function which is integrated to obtain the second stage utiltiy.
```julia-repl
julia> Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
```

* ```f1(Con,Stat, t::Float64, Para::Dict)``` and ```f2(Con,Stat, t::Float64, s::Float64, Para::Dict)``` are the function of the right hand side of the state dynamics.
```julia-repl
julia> f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1],-Para["eta"]*Stat[2]]
julia> f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1], 0]
```

* ```g(Con,Stat,t::Float64, Para::Dict)``` describes the transition function between the state variables in the two stages.
```julia-repl
julia> g(Con,Stat,t::Float64, Para::Dict) = [(1-Para["gamma"])*Stat[1],Para["eta"]*Stat[2]]
```

* ```S1(Stat, Para::Dict)``` and ```S2(Stat_dist, Para::Dict)``` are the salvage-function in the first and second stage. Since we have no salvage value in our toy-example, we could omit them as entries in ```TwoStageOptimisation```. However, we can also specifiy and enter them.
```julia-repl
julia> S1(Stat, Para::Dict) = 0
julia> S2(Stat_dist, Para::Dict) = 0
```

## Solve the model

Finally we can solve the model. But first, we specify, that we want to give the initial guess for the control profiles to the algorithm
```julia-repl
julia> MyPara["LoadInits"] = true
```

Now we define the profiles. We can pick the dimensions across time how big we like, the controls get interpolated in the solution process. Here we decided for 100 points. We define the control profiles as entries in the ```Results``` dictionary.

```julia-repl
julia> Results = Dict()
julia> Results["Con"] = 4.0*ones(1,100,1)
julia> Results["Con_dist"] = 4.0*ones(100,100,1)
julia> Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand1 = U, 
                                AggregationIntegrand2 = Q,
                                StateDynamic1 = f1,
                                StateDynamic2 = f2, 
                                Shock12 = g,
                                SalvageFunction1=S1,
                                SalvageFunction2=S2)
```

Finally we can save the results in location of choice and also plot the final results. We have the option to save the plots, and they are automatically save

```julia-repl
julia> SaveResults(Results,"CapitalAccumulation")
julia> PlotResults(Results,SavePlot = true)
```