using TwoStageOptimalControl
using Test


@testset "TestFirstStageSolution" begin
    # Test code for a test model without shock
    ResultsFirstStage = LoadResults("test/TestModelFirstStage")
    MyPara = deepcopy(ResultsFirstStage["Para"])
    MyPara["hstep"] = MyPara["hstepStart"]
    MyPara["PlotResultsStart"] = false
    MyPara["PlotResultsIntermediate"] = false
    MyPara["PlotResultsFinal"] = false
    MyPara["SearchDISP"] = false
    delete!(MyPara,["nTime","nAge","tmesh","amesh"])

    U(Con, Stat, t::Float64, Para::Dict) = Stat[1]^0.5 - Con[1]^2
    Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = 0
    f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1]]
    f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = []
    g(Con,Stat,t::Float64, Para::Dict) = []
    S1(Stat, Para::Dict) = 0
    S2(Stat_dist, Para::Dict) = 0
    
    Results = Dict()
    Results["Con"] = 3.0*ones(1,100,1)
    Results["Con_dist"] = 3.0*ones(1,2,0)
    Results = TwoStageOptimalControl.TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                    ObjectiveIntegrand1 = U, 
                                    AggregationIntegrand2 = Q,
                                    StateDynamic1 = f1,
                                    StateDynamic2 = f2, 
                                    Shock12 = g,
                                    SalvageFunction1 = S1,
                                    SalvageFunction2 = S2,
    )

    maxControlDifference = maximum(abs.(Results["Con"]-ResultsFirstStage["Con"]))
    maxStateDifference = maximum(abs.(Results["Stat"]-ResultsFirstStage["Stat"]))
    maxCostateDifference = maximum(abs.(Results["CoStat"]-ResultsFirstStage["CoStat"]))

    @test maxControlDifference < 1e-5
    @test maxStateDifference < 1e-5
    @test maxCostateDifference < 1e-5
end

@testset "TestModel" begin
    # Test code for a test model without shock
    ResultsBenchmark = LoadResults("test/TestModelBenchmark")
    MyPara = deepcopy(ResultsBenchmark["Para"])
    MyPara["hstep"] = MyPara["hstepStart"]
    MyPara["PlotResultsStart"] = false
    MyPara["PlotResultsIntermediate"] = false
    MyPara["PlotResultsFinal"] = false
    MyPara["SearchDISP"] = false
    delete!(MyPara,["nTime","nAge","tmesh","amesh"])

    U(Con, Stat, t::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
    Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
    f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1],
                                                    -0.5*Stat[2]]
    f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1], 
                                                     0]
    g(Con,Stat,t::Float64, Para::Dict) = [0.8*Stat[1],0.5*Stat[2]]
    S1(Stat, Para::Dict) = Stat[2]*0
    S2(Stat_dist, Para::Dict) = Stat_dist[2]*0

    Results = Dict()
    Results["Con"] = 4.0*ones(1,100,1)
    Results["Con_dist"] = 4.0*ones(100,100,1)
    Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand1 = U, 
                                AggregationIntegrand2 = Q,
                                StateDynamic1 = f1,
                                StateDynamic2 = f2, 
                                Shock12 = g)

    maxControlDifference = max(maximum(abs.(Results["Con"]-ResultsBenchmark["Con"])),
                               maximum(abs.(Results["Con_dist"]-ResultsBenchmark["Con_dist"])))
    maxStateDifference = max(maximum(abs.(Results["Stat"]-ResultsBenchmark["Stat"])),
                               maximum(abs.(Results["Stat_dist"]-ResultsBenchmark["Stat_dist"])))
    maxCostateDifference = max(maximum(abs.(Results["CoStat"]-ResultsBenchmark["CoStat"])),
                               maximum(abs.(Results["CoStat_dist"]-ResultsBenchmark["CoStat_dist"])))

    @test maxControlDifference < 1e-5
    @test maxStateDifference < 1e-5
    @test maxCostateDifference < 1e-5
end

@testset "TestModel 2" begin
    # Test code for a capital accumulation model without shock

    ResultsBenchmark = LoadResults("test/TestModelBenchmark2")

    MyPara = deepcopy(ResultsBenchmark["Para"])
    MyPara["hstep"] = MyPara["hstepStart"]
    MyPara["PlotResultsStart"] = false
    MyPara["PlotResultsIntermediate"] = false
    MyPara["PlotResultsFinal"] = false
    MyPara["SearchDISP"] = false
    delete!(MyPara,["nTime","nAge","tmesh","amesh"])

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

    Results = Dict()
    Results["Con"] = 0.0*ones(1,100,1)
    Results["Con_dist"] = 0.0*ones(100,100,1)
    
    Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand1 = U, 
                                AggregationIntegrand2 = Q,
                                StateDynamic1 = f1,
                                StateDynamic2 = f2, 
                                Shock12 = g,
                                SalvageFunction1=S1,
                                SalvageFunction2=S2)

    maxControlDifference = max(maximum(abs.(Results["Con"]-ResultsBenchmark["Con"])),
                               maximum(abs.(Results["Con_dist"]-ResultsBenchmark["Con_dist"])))
    maxStateDifference = max(maximum(abs.(Results["Stat"]-ResultsBenchmark["Stat"])),
                               maximum(abs.(Results["Stat_dist"]-ResultsBenchmark["Stat_dist"])))
    maxCostateDifference = max(maximum(abs.(Results["CoStat"]-ResultsBenchmark["CoStat"])),
                               maximum(abs.(Results["CoStat_dist"]-ResultsBenchmark["CoStat_dist"])))

    @test maxControlDifference < 1e-5
    @test maxStateDifference < 1e-5
    @test maxCostateDifference < 1e-5
end
