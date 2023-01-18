using TwoStageOptControl
using Test

@testset "Capital accumulation" begin
    # Test code for a capital accumulation model without shock

    ResultsBenchmark = LoadResults("test/CapitalAccumulationBenchark")

    MyPara = copy(ResultsBenchmark["Para"])
    U(Con, Stat, t::Float64, Para::Dict) = Stat[2] * ((Para["A"] - Stat[1])*Stat[1] - Para["b"]*Con[1] - Para["c"]/2*Con[1]^2)
    Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = Stat[2] * ((Para["A"] - Stat[1])*Stat[1] - Para["b"]*Con[1] - Para["c"]/2*Con[1]^2)
    f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1], -Para["eta"]*Stat[2]]
    f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - Para["delta"]*Stat[1], 0]
    g(Con,Stat,t::Float64, Para::Dict) = [Stat[1],Para["eta"]*Stat[2]]
    S1(Stat, Para::Dict) = 0
    S2(Stat_dist, Para::Dict) = 0
    @test my_f(2,1) == 7

end


@testset "Test model" begin
    # Test code for a test model without shock

    ResultsBenchmark = LoadResults("test/TestModelBenchark")

    MyPara = copy(ResultsBenchmark["Para"])
    U(Con, Stat, t::Float64, Para::Dict) = Stat[2]*(Stat[1]^0.5 - Con[1]^2)
    Q(Con,Stat, t::Float64,s::Float64, Para::Dict) = 1*Stat[2]*(Stat[1]^0.5 - Con[1]^2)
    f1(Con,Stat, t::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1],-0.5*Stat[2]]
    f2(Con,Stat, t::Float64, s::Float64, Para::Dict) = [Con[1] - 0.1*Stat[1], 0]
    g(Con,Stat,t::Float64, Para::Dict) = [Stat[1],0.5*Stat[2]]
    S1(Stat, Para::Dict) = 0
    S2(Stat_dist, Para::Dict) = 0
    @test my_f(2,1) == 7

end