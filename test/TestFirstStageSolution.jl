#using TwoStageOptimalControl

U(Con, Stat, t::Float64, Para::Dict) = Stat[1]^0.5 - Con[1]^2
Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = 0

f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1]]
f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = []

g(Con,Stat,t::Float64, Para::Dict) = []

S1(Stat, Para::Dict) = 0
S2(Stat_dist, Para::Dict) = 0


MyPara = Dict()
MyPara["T"] = 1
MyPara["hstep"] = 0.05

MyPara["nCon"] = 1
MyPara["nCon_dist"] = 0
MyPara["nStat"] = 1
MyPara["nStat_dist"] = 0
MyPara["InitStat"] = [0.1]
MyPara["ConMin"] = (t) -> [0.0]

MyPara["InitLineStep"] = 1e-6
MyPara["UpperLineStep"] = 1e-2
MyPara["hLowBound"] = 0.007
MyPara["PlotResultsIntermediateFrequency"] = 150
MyPara["OptiType"] = "Newton-Raphson"
MyPara["PlotResultsWaitForKey"] = true

MyPara["LoadInits"] = true
#MyPara["FixedControl1"] = [1]

Results = Dict()
Results["Con"] = 4.0*ones(1,100,1)
Results["Con_dist"] = 4.0*ones(1,2,0)
Results = Main.TwoStageOptimalControl.TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand1 = U, 
                                AggregationIntegrand2 = Q,
                                StateDynamic1 = f1,
                                StateDynamic2 = f2, 
                                Shock12 = g,
                                SalvageFunction1 = S1,
                                SalvageFunction2 = S2,
)

SaveResults(Results,"test/TestModelBenchmark")


