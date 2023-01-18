U(Con, Stat, t::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = 1*Stat[2]*(Stat[1]^0.5 - Con[1]^2)

f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1],
                                        -0.5*Stat[2]]
f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1],
                                                    0]

g(Con,Stat,t::Float64, Para::Dict) = [Stat[1],0.5*Stat[2]]

S1(Stat, Para::Dict) = 0
S2(Stat_dist, Para::Dict) = 0


MyPara = Dict()
MyPara["T"] = 1
MyPara["hstep"] = 0.05

MyPara["nCon"] = 1
MyPara["nCon_dist"] = 1
MyPara["nStat"] = 2
MyPara["nStat_dist"] = 2
MyPara["InitStat"] = [0.1,1.0]
MyPara["ConMin"] = [0.0]
MyPara["Con_distMin"] = [0.0]

MyPara["InitLineStep"] = 1e-5
MyPara["UpperLineStep"] = 1e-1
MyPara["hLowBound"] = 0.01
MyPara["PlotResultsIntermediateFrequency"] = 30

MyPara["LoadInits"] = true

Results = Dict()
Results["Con"] = 4.0*ones(1,20,1)
Results["Con_dist"] = 4.0*ones(20,20,1)
Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand2 = U, 
                                AggregationFunction2 = Q,
                                StateDynamic_1_2 = f1,
                                StateDynamic_2_2 = f2, 
                                Shock2 = g,
                                SalvageFunction_1_2=S1,
                                SalvageFunction_2_2=S2)

SaveResults(Results,"test/TestModelBenchmark")


