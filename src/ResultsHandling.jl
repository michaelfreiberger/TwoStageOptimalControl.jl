
"""
    AssignResults(Results::Dict, Con::Array{Float64,3}, Stat::Array{Float64,3}, Con_dist::Array{Float64,3}, Stat_dist::Array{Float64,3}, Stat_agg::Array{Float64,3}, CoStat::Array{Float64,3}, CoStat_dist::Array{Float64,3}, dHam::Array{Float64,3}, dHam_dist::Array{Float64,3}, Para::Dict)

    Assign the control and state variables to the Results dictionary
"""
function AssignResults(Results::Dict, Con::Array{Float64,3}, Stat::Array{Float64,3}, Con_dist::Array{Float64,3}, Stat_dist::Array{Float64,3}, Stat_agg::Array{Float64,3}, CoStat::Array{Float64,3}, CoStat_dist::Array{Float64,3}, dHam::Array{Float64,3}, dHam_dist::Array{Float64,3}, Para::Dict)

    Results["Con"] = Con
    Results["Stat"] = Stat
    Results["CoStat"] = CoStat
    Results["Stat_agg"] = Stat_agg
    Results["dHam"] = dHam
    Results["Con_dist"] = Con_dist
    Results["Stat_dist"] = Stat_dist
    Results["CoStat_dist"] = CoStat_dist
    Results["dHam_dist"] = dHam_dist
    Results["Para"] = Para
    Results["ObjectiveValue"] = ObjectValue(Con,Stat,Con_dist,Stat_dist,Stat_agg,Para)

    return
end

"""
    PlotResults(Results2::Dict;SavePlot=false,Display=true,sizeX=600,sizeY=400)

    Basic plots of all controls + gradient and state + costate variables
"""
function PlotResults(Results2::Dict;SavePlot=false,Display=true,sizeX=600,sizeY=400)

    Results = deepcopy(Results2)
    Con_dist = Results["Con_dist"]
    Con = Results["Con"]
    Stat_dist = Results["Stat_dist"]
    Stat = Results["Stat"]
    Stat_agg = Results["Stat_agg"]
    CoStat_dist = Results["CoStat_dist"]
    CoStat = Results["CoStat"]
    dHam_dist = Results["dHam_dist"]
    dHam = Results["dHam"]
    Para = Results["Para"]

    PLOTS = Dict();

    nPlotHelp = min(Para["nTime"],Para["nVintagePlot"])
    ylab = string("")

    if SavePlot == true
        pyplot();
        ControlLabels = Para["ControlLabelsLatex"]
        StateLabels = Para["StateLabelsLatex"]
        CoStateLabels = Para["CoStateLabelsLatex"]
        ControlDistLabels = Para["ControlDistLabelsLatex"]
        StateDistLabels = Para["StateDistLabelsLatex"]
        CoStateDistLabels = Para["CoStateDistLabelsLatex"]
    else
        gr();
        ControlLabels = Para["ControlLabelsSimple"]
        StateLabels = Para["StateLabelsSimple"]
        CoStateLabels = Para["CoStateLabelsSimple"]
        ControlDistLabels = Para["ControlDistLabelsSimple"]
        StateDistLabels = Para["StateDistLabelsSimple"]
        CoStateDistLabels = Para["CoStateDistLabelsSimple"]
    end

	colors = get(colorschemes[:rainbow_bgyr_35_85_c72_n256],range(0,stop=1,length=nPlotHelp))

    #------------------------------------------------------------------------------------
    #   Plot of non-distributed Controls
    PLOTS["Controls"] = plot(Para["tmesh"],Con[1,:,:],lw=2,title = "Controls",legend=:best,label=reshape(ControlLabels,1,:),size = (sizeX, sizeY))

    #   Plot of Gradient w.r.t. non-distributed Controls
    PLOTS["Gradients"] = plot(Para["tmesh"],dHam[1,:,:],lw=2,title = "Gradients",legend=:best,label=reshape(ControlLabels,1,:),size = (sizeX, sizeY))

    #   Plot of non-distributed States
    PLOTS["States"] = plot(Para["tmesh"],Stat[1,:,:],lw=2,title = "States",legend=:best,label=reshape(StateLabels,1,:),size = (sizeX, sizeY))

    #   Plot of non-distributed Co-states
    PLOTS["CoStates"] = plot(Para["tmesh"],CoStat[1,:,:],lw=2,title = "CoStates",legend = :best,label=reshape(CoStateLabels,1,:),size = (sizeX, sizeY))

    #------------------------------------------------------------------------------------
    #   Plot of distributed controls
    for k = 1:Para["nCon_dist"]
        PLOTS[string("Con_dist",k)] = plot(title = string("Distributed Control ",ControlDistLabels[k]),size = (sizeX, sizeY),ylabel = ylab)
        for ii=1:nPlotHelp
            index = round(Int,(ii-1)/(nPlotHelp-1)*(Para["nTime"]-1) + 1)
            PLOTS[string("Con_dist",k)] = plot!(Para["tmesh"][index:end],Con_dist[index,index:end,k],lw = 2,color=colors[ii],label=false)
        end
        PLOTS[string("Con_dist",k)] = plot!(Para["tmesh"],diag(Con_dist[:,:,k]),lw = 2,color=:black,line=(:dash,2),
                                            label=string("Initial second stage values ",ControlDistLabels[k],"(t,t)"))
        PLOTS[string("Con_dist",k)] = plot!(Para["tmesh"],Con_dist[1,:,k],lw = 2,color=:grey,line=(:dot,2),
                                            label=string("First stage values ",ControlDistLabels[k],"(t)"))
    end

    #------------------------------------------------------------------------------------
    #   Plot of distributed gradient
    for k = 1:Para["nCon_dist"]
        PLOTS[string("Grad_dist",k)] = plot(title = string("Gradient of distributed Control ",ControlDistLabels[k]),size = (sizeX, sizeY))
        for ii=1:nPlotHelp
            index = round(Int,(ii-1)/Para["hstep"]*Para["T"]/(nPlotHelp-1) + 1)
            PLOTS[string("Grad_dist",k)] = plot!(Para["tmesh"][index:end],dHam_dist[index,index:end,k],lw = 2,color=colors[ii],label=false)
        end
        PLOTS[string("Grad_dist",k)] = plot!(Para["tmesh"],dHam_dist[1,:,k],lw = 2,color=:black,line=(:dot,2),
                                            label=string("First stage values ",ControlDistLabels[k],"(t)"))
    end

    #------------------------------------------------------------------------------------
    #   Plot of distributed states
    for k = 1:Para["nStat_dist"]
        PLOTS[string("Stat_dist",k)] = plot(title = string("Distributed State ",StateDistLabels[k]),size = (sizeX, sizeY))
        for ii=1:nPlotHelp
            index = round(Int,(ii-1)/Para["hstep"]*Para["T"]/(nPlotHelp-1) + 1)
            PLOTS[string("Stat_dist",k)] = plot!(Para["tmesh"][index:end],Stat_dist[index,index:end,k],lw = 2,color=colors[ii],label=false)
        end
        PLOTS[string("Stat_dist",k)] = plot!(Para["tmesh"],diag(Stat_dist[:,:,k]),lw = 2,color=:black,line=(:dash,2),
                                            label=string("Initial second stage values ",StateDistLabels[k],"(t,t)"))
    end

    #------------------------------------------------------------------------------------
    #   Plot of distributed costates
    for k = 1:Para["nStat_dist"]
        PLOTS[string("CoStat_dist",k)] = plot(title = string("Distributed CoState ",CoStateDistLabels[k]),size = (sizeX, sizeY))
        for ii=1:nPlotHelp
            index = round(Int,(ii-1)/Para["hstep"]*Para["T"]/(nPlotHelp-1) + 1)
            PLOTS[string("CoStat_dist",k)] = plot!(Para["tmesh"][index:end],CoStat_dist[index,index:end,k],lw = 2,color=colors[ii],label=false)
        end
        PLOTS[string("CoStat_dist",k)] = plot!(Para["tmesh"],diag(CoStat_dist[:,:,k]),lw = 2,color=:black,line=(:dash,2),
                                            label=string("Initial second stage values ",StateDistLabels[k],"(t,t)"))
        PLOTS[string("CoStat_dist",k)] = plot!(Para["tmesh"],CoStat_dist[1,:,k],lw = 2,color=:grey,line=(:dot,2),
                                            label=string("First stage values ",StateDistLabels[k],"(t)"))                                        
    end


    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    if SavePlot == true
        mkpath(string(Para["SavePlotPath"],"/BasePlots"))
        for str in sort(collect(keys(PLOTS)))
            savefig(PLOTS[str],string(Para["SavePlotPath"],"/BasePlots/",str,".pdf"))
        end
    end

    if Display == true
        l = @layout [a;b]
        display(plot(PLOTS["Controls"],PLOTS["Gradients"],layout = l,size=[600,800]))
        l = @layout [a;b]
        display(plot(PLOTS["States"],PLOTS["CoStates"],layout = l,size=[600,800]))

        for k=1:Para["nCon_dist"]
            l = @layout [a;b]
            display(plot(PLOTS[string("Con_dist",k)],PLOTS[string("Grad_dist",k)],layout = l,size=[600,800]))
        end
        for k=1:Para["nStat_dist"]
           l = @layout [a;b]
           display(plot(PLOTS[string("Stat_dist",k)],PLOTS[string("CoStat_dist",k)],layout = l,size=[600,800]))
        end
    end
    return PLOTS
end


"""
    SaveResults(Results::Dict,filepath)

    Save the Results Dictionary to filepath
"""
function SaveResults(Results::Dict,filepath)

    Results["Para"]["SavePlotPath"] = filepath

    open(string(filepath,".out"),"w") do f
        serialize(f,Results)
    end
end

"""
    LoadResults(filepath)

    Load the Results Dictionary from the filepath
"""
function LoadResults(filepath)
    Results = Dict()
    open(string(filepath,".out")) do f
        Results = deserialize(f)
    end
    return Results
end