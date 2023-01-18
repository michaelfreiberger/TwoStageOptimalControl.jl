
U(Con, Stat, t::Float64, Para::Dict) = Stat[2] * ((Para["A"] - Stat[1])*Stat[1] - Para["b"]*Con[1] - Para["c"]/2*Con[1]^2)
Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = Stat[2] * ((Para["A"] - Stat[1])*Stat[1] - Para["b"]*Con[1] - Para["c"]/2*Con[1]^2)

f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1],
                                        -Para["eta"]*Stat[2]]
f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1],
                                                    0]

g(Con,Stat,t::Float64, Para::Dict) = [Stat[1],
                                      Para["eta"]*Stat[2]]

S1(Stat, Para::Dict) = 0
S2(Stat_dist, Para::Dict) = 0

MyPara = Dict()
MyPara["A"] = 60
MyPara["b"] = 100
MyPara["c"] = 10
MyPara["delta"] = 0.1
MyPara["rho"] = 0.1
MyPara["eta"] = 0.05

MyPara["T"] = 10
MyPara["hstep"] = 0.5
MyPara["nCon"] = 1
MyPara["nCon_dist"] = 1
MyPara["nStat"] = 2
MyPara["nStat_dist"] = 2
MyPara["InitStat"] = [0.1,1.0]
MyPara["ConMin"] = [0.0]
MyPara["Con_distMin"] = [0.0]


MyPara["OptiType"] = "Newton-Raphson"
MyPara["ProbIndex"] = 2

MyPara["InitLineStep"] = 1e-5
MyPara["UpperLineStep"] = 1e-2
MyPara["hLowBound"] = 0.1
MyPara["PlotResultsIntermediateFrequency"] = 150

MyPara["LoadInits"] = true

Results = Dict()
Results["Con"] = 6.0*ones(1,10,1)
Results["Con_dist"] = 6.0*ones(10,10,1)
Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand2 = U, 
                                AggregationFunction2 = Q,
                                StateDynamic_1_2 = f1,
                                StateDynamic_2_2 = f2, 
                                Shock2 = g,
                                SalvageFunction_1_2=S1,
                                SalvageFunction_2_2=S2)

SaveResults(Results,"test/CapitalAccumulationBenchmark")
