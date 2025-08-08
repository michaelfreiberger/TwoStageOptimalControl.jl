using TwoStageOptimalControl

U(Con, Stat, t::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
Q(Con, Stat, t::Float64,s::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)

f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1],
                                        -0.5*Stat[2]]
f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1],
                                                    0]

g(Con,Stat,t::Float64, Para::Dict) = [0.8*Stat[1],0.5*Stat[2]]

S1(Stat, Para::Dict) = Stat[2]*0
S2(Stat_dist, Para::Dict) = Stat_dist[2]*0


MyPara = Dict()
MyPara["T"] = 1
MyPara["hstep"] = 0.05

MyPara["nCon"] = 1
MyPara["nCon_dist"] = 1
MyPara["nStat"] = 2
MyPara["nStat_dist"] = 2
MyPara["InitStat"] = [0.1,1.0]
MyPara["ConMin"] = (t) -> [0.0]
MyPara["Con_distMin"] = (t) -> [0.0]

MyPara["InitLineStep"] = 1e-6
MyPara["UpperLineStep"] = 1e-2
MyPara["hLowBound"] = 0.007
MyPara["PlotResultsIntermediateFrequency"] = 150
MyPara["OptiType"] = "Newton-Raphson"
MyPara["PlotResultsWaitForKey"] = false

MyPara["LoadInits"] = true

Results = Dict()
Results["Con"] = 4.0*ones(1,100,1)
Results["Con_dist"] = 4.0*ones(100,100,1)
Results = TwoStageOptimalControl.TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand1 = U, 
                                AggregationIntegrand2 = Q,
                                StateDynamic1 = f1,
                                StateDynamic2 = f2, 
                                Shock12 = g,
                                SalvageFunction1 = S1,
                                SalvageFunction2 = S2,
)

TwoStageOptimalControl.SaveResults(Results,"test/TestModelBenchmark")