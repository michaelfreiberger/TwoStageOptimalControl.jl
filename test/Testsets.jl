using TwoStageOptimalControl
using Test


@testset "TestModel" begin
    # Test code for a test model without shock
    ResultsBenchmark = LoadResults("TestModelBenchmark")
    MyPara = deepcopy(ResultsBenchmark["Para"])
    MyPara["hstep"] = 0.05
    MyPara["OptiType"] = "Newton-Raphson"
    MyPara["PlotResultsStart"] = false
    MyPara["PlotResultsIntermediate"] = false
    MyPara["PlotResultsFinal"] = false
    MyPara["SearchDISP"] = false
    delete!(MyPara,["nTime","nAge","tmesh","amesh"])

    U(Con, Stat, t::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
    Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = 1*Stat[2]*(Stat[1]^0.5 - Con[1]^2)
    f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1],-0.5*Stat[2]]
    f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1], 0]
    g(Con,Stat,t::Float64, Para::Dict) = [Stat[1],0.5*Stat[2]]
    
    Results = Dict()
    Results["Con"] = 3.0*ones(1,100,1)
    Results["Con_dist"] = 3.0*ones(100,100,1)
    
    Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand2 = U, 
                                AggregationFunction2 = Q,
                                StateDynamic_1_2 = f1,
                                StateDynamic_2_2 = f2, 
                                Shock2 = g)

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


@testset "Capital accumulation" begin
    # Test code for a capital accumulation model without shock

    ResultsBenchmark = LoadResults("CapitalAccumulationBenchmark")

    MyPara = deepcopy(ResultsBenchmark["Para"])
    MyPara["hstep"] = 0.5
    MyPara["OptiType"] = "Newton-Raphson"
    MyPara["PlotResultsStart"] = false
    MyPara["PlotResultsIntermediate"] = false
    MyPara["PlotResultsFinal"] = false
    MyPara["SearchDISP"] = false
    delete!(MyPara,["nTime","nAge","tmesh","amesh"])

    U(Con, Stat, t::Float64, Para::Dict) = Stat[2] * ((Para["A"] - Stat[1])*Stat[1] - Para["b"]*Con[1] - Para["c"]/2*Con[1]^2)
    Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = Stat[2] * ((Para["A"] - Stat[1])*Stat[1] - Para["b"]*Con[1] - Para["c"]/2*Con[1]^2)
    f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1], -Para["eta"]*Stat[2]]
    f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1], 0]
    g(Con,Stat,t::Float64, Para::Dict) = [Stat[1],Para["eta"]*Stat[2]]
    
    Results = Dict()
    Results["Con"] = 6.0*ones(1,100,1)
    Results["Con_dist"] = 6.0*ones(100,100,1)
    
    Results = TwoStageOptimisation(Results = Results,UserParameters = MyPara,
                                ObjectiveIntegrand2 = U, 
                                AggregationFunction2 = Q,
                                StateDynamic_1_2 = f1,
                                StateDynamic_2_2 = f2, 
                                Shock2 = g)

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
