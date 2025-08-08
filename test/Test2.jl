using TwoStageOptimalControl

U(Con, Stat, t::Float64, Para::Dict) = Stat[2] * ((Para["A1"] - Stat[1])*Stat[1] - Para["b1"]*Con[1] - Para["c1"]/2*Con[1]^2)
Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = Stat[2] * ((Para["A2"] - Stat[1])*Stat[1] - Para["b2"]*Con[1] - Para["c2"]/2*Con[1]^2)

f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1],
                                        -Para["eta"]*Stat[2]]
f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1],
                                                    0]

g(Con,Stat,t::Float64, Para::Dict) = [Stat[1],
                                      Para["eta"]*Stat[2]]

S1(Stat, Para::Dict) = 0
S2(Stat_dist, Para::Dict) = 0

MyPara = Dict()
MyPara["A1"] = 60
MyPara["A2"] = 80
MyPara["b1"] = 100
MyPara["b2"] = 120
MyPara["c1"] = 10
MyPara["c2"] = 15
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
MyPara["ConMin"] = (t) -> [0.0]
MyPara["Con_distMin"] = (t) -> [0.0]


MyPara["OptiType"] = "Newton-Raphson"
MyPara["ProbIndex"] = 2
MyPara["InitLineStep"] = 1e-5
MyPara["UpperLineStep"] = 1e-1
MyPara["hLowBound"] = 0.05
MyPara["PlotResultsIntermediateFrequency"] = 500
MyPara["OptiType"] = "Newton-Raphson"

MyPara["LoadInits"] = true

Results = Dict()
Results["Con"] = 0.0*ones(1,100,1)
Results["Con_dist"] = 0.0*ones(100,100,1)
Results = TwoStageOptimalControl.TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand1 = U, 
                                AggregationIntegrand2 = Q,
                                StateDynamic1 = f1,
                                StateDynamic2 = f2, 
                                Shock12 = g,
                                SalvageFunction1=S1,
                                SalvageFunction2=S2)

TwoStageOptimalControl.SaveResults(Results,"test/TestModelBenchmark2")
