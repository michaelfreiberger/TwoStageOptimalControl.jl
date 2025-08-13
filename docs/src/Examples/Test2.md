# Second optimisation problem

Here we investigate a second simple capital accumulation problem. A firm optimises profits which arise from a quadratic production function and a quadratic cost function. After a random shock, the production and cost parameters change.

```math
\begin{aligned}
\max_{I_1(t),I_2(t,s)} &\int_0^T e^{\rho t}\cdot\left\{ Z_1(t)\cdot \Big[ (A  - K_1(t))*K_1(t) - b_1*I_1(t) - c_1/2*I_1(t)^2\Big] + Q(t) \right\} dt \\
s.t.\hspace{1cm} \dot{K_1}(t) &= I_1(t) - \delta K_1(t) &&,\qquad K_1(0) = K_1^0 \\
\dot{Z_1}(t) &= - \eta\cdot Z_1(t)                      &&,\qquad Z_1(0) = 1 \\
\dot{K_2}(t,s) &= I_2(t,s) - \delta K_2(t,s)            &&,\qquad K_2(s,s) = (1-\gamma)\cdot K_1(s) \\
\dot{Z_2}(t,s) &= 0                                     &&,\qquad Z_2(s,s) = \eta\cdot Z_1(s) \\
Q(t) &= \int_0^t Z_2(t,s)\cdot \Big[ (A_2  - K_2(t,s))*K_2(t,s) - b_2*I_2(t,s) - c_2/2*I_2(t,s)^2\Big] ds
\end{aligned}
```

```julia-repl
julia> using TwoStageOptimalControl

julia> U(Con, Stat, t::Float64, Para::Dict) = Stat[2] * ((Para["A1"] - Stat[1])*Stat[1] - Para["b1"]*Con[1] - Para["c1"]/2*Con[1]^2)
julia> Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = Stat[2] * ((Para["A2"] - Stat[1])*Stat[1] - Para["b2"]*Con[1] - Para["c2"]/2*Con[1]^2)
julia> f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1],   -Para["eta"]*Stat[2]]
julia> f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1],0]
julia> g(Con,Stat,t::Float64, Para::Dict) = [Stat[1],Para["eta"]*Stat[2]]
julia> S1(Stat, Para::Dict) = 0
julia> S2(Stat_dist, Para::Dict) = 0
```

```julia-repl
julia> MyPara = Dict()
julia> MyPara["A1"] = 60
julia> MyPara["A2"] = 80
julia> MyPara["b1"] = 100
julia> MyPara["b2"] = 120
julia> MyPara["c1"] = 10
julia> MyPara["c2"] = 15
julia> MyPara["delta"] = 0.1
julia> MyPara["rho"] = 0.1
julia> MyPara["eta"] = 0.05

julia> MyPara["T"] = 10
julia> MyPara["hstep"] = 0.5
julia> MyPara["nCon"] = 1
julia> MyPara["nCon_dist"] = 1
julia> MyPara["nStat"] = 2
julia> MyPara["nStat_dist"] = 2
julia> MyPara["InitStat"] = [0.1,1.0]
julia> MyPara["ConMin"] = t->[0.0]
julia> MyPara["Con_distMin"] = t->[0.0]

julia> MyPara["OptiType"] = "Newton-Raphson"
julia> MyPara["ProbIndex"] = 2

julia> MyPara["InitLineStep"] = 1e-5
julia> MyPara["UpperLineStep"] = 1e-2
julia> MyPara["hLowBound"] = 0.1
julia> MyPara["PlotResultsIntermediateFrequency"] = 150

julia> MyPara["LoadInits"] = true
```

```julia-repl
julia> Results = Dict()
julia> Results["Con"] = 6.0*ones(1,10,1)
julia> Results["Con_dist"] = 6.0*ones(10,10,1)
julia> Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand = U, 
                                AggregationFunction = Q,
                                StateDynamic_1 = f1,
                                StateDynamic_2 = f2, 
                                Shock = g,
                                SalvageFunction_1=S1,
                                SalvageFunction_2=S2)

julia> SaveResults(Results,"test/CapitalAccumulationBenchmark")
```
