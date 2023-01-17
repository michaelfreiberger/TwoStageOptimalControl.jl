
"""
    ParametersBasics(Para,UserPara)

    Initialize the basic parameters necessary for the model
"""
function ParametersBasics(Para,UserPara)
    #= --------------------------------------------------------------
        Parameters for structure/dimensions of the model
    -------------------------------------------------------------- =#
    
    Para2 = Dict()

    # time-horizont
    Para2["T"] = 1.0                       
    
    # stepsize
    Para2["hstep"] = 0.05                     

    # discount rate
    Para2["rho"] = 0.0
    
    # number of controls in first stage
    Para2["nCon"] = 2                   
    # number of controls in second stage          
    Para2["nCon_dist"] = 2          
    # number of state varibales in first stage
    Para2["nStat"] = 3           
    # number of state varibales in second stage                       
    Para2["nStat_dist"] = 3 
    
    # index of state variable describing the auxiliary variables in both stages
    Para2["ProbIndex"] = 1

    # Overwrite values if given by the User
    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end

"""
    ParametersGrids(Para,UserPara)

    Initialize the parameters associated with the time grids
"""
function ParametersGrids(Para,UserPara)
    Para2 = Dict()
    # number of gridpoints
    Para2["nTime"] = floor(Int,Para["T"]/Para["hstep"]) + 1 
    Para2["nAge"] = floor(Int,Para["T"]/Para["hstep"]) + 1 

    # grid over the time horizon
    Para2["tmesh"] = collect(0:Para["hstep"]:Para["T"])     
    Para2["amesh"] = collect(0:Para["hstep"]:Para["T"]) 

    # Overwrite values if given by the User
    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end

"""
    ParametersControls(Para,UserPara)

    Initialize the parameters associated with the control variables
"""
function ParametersControls(Para,UserPara)
    Para2 = Dict()
    # Loading specific initial guesses for the control variables
    Para2["LoadInits"] = false                   
    # Indices of controls, which should be disabled (Stage 1)
    Para2["Reduce1"] = []                       
    # Indices of controls, which should be disabled (Stage 2)
    Para2["Reduce2"] = []       
    # Lower and Upper bounds for the control variables
    Para2["ConMin"] = -Inf*ones(Para["nCon"])
    Para2["ConMax"] =  Inf*ones(Para["nCon"])
    Para2["Con_distMin"] = -Inf*ones(Para["nCon_dist"])
    Para2["Con_distMax"] = Inf*ones(Para["nCon_dist"])

    # Overwrite values if given by the User
    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end

"""
    ParametersAlgorithm(Para,UserPara)

    Initialize the parameters for the optimisation algorithm
"""
function ParametersAlgorithm(Para,UserPara)
    #= ------------------------------------------------------------
        Parameters for the optimisation algorithm
    ------------------------------------------------------------ =#
    Para2 = Dict()    
    # Adjustment of gradient of Hamiltonian
    #   - "Gradient" uses the standard gradient of the hamiltonian 
    #   - "ProbAdjust" scales the gradient by the inverse of the auxiliary variables
    #   - "Newton-Raphson" scales the gradient by the inverse of the Hessian
    Para2["OptiType"] = "Gradient" 

    Para2["IntegrationOrder"] = 4
    Para2["ParallelComputing"] = false
    
    # Lower bound for the step-size
    Para2["hLowBound"] = Para["hstep"]*0.9 

    # Initial stepsize for linesearch
    Para2["InitLineStep"] = 1e-6          
    # Factor for increase in the linesearch stepsize
    Para2["LineSearchStepIncrease"] = 0.25
    # Upper bound for linesearch stepsize
    Para2["UpperLineStep"] = 1e1           
    # Lower bound for step in gradient search -> termination if undercut
    Para2["stepLowBound"] = 1e-8   
    # Lower bound for the gradient -> termination if undercut
    Para2["GradBound"] = 1e-5
    # Max iterations with same step-size h
    Para2["MaxOptiIter"] = 10000                   
    # Max iterations along the same gradient 
    Para2["MaxLineIter"] = 30                      

    # Indicator if gradients are smoothed
    Para2["GradSmooth"] = false              
    # Indicator if controls are smoothed  
    Para2["ConSmooth"] = false               

    # Indicator if intermediate results shall be displayed
    Para2["SearchDISP"] = true                 

    # Iteration counters
    Para2["LineIter"] = 0                      
    Para2["OptiIter"] = 0
    Para2["HstepIter"] = 0

    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end

function ParametersInitStates(Para,UserPara)
    #= ------------------------------------------------------------
        Iinitial and boundary data
    ------------------------------------------------------------ =#
    Para2 = Dict()
    Para2["InitStat"] = zeros(1,1,Para["nStat"])
    Para2["InitStat_dist"] = zeros(Para["nAge"],1,Para["nStat_dist"])
    Para2["BoundStat_dist"] = zeros(1,Para["nTime"],Para["nStat_dist"])

    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end

"""
    ParametersPlots(Para,UserPara)

    Initialize parameters relevant for plotting
"""
function ParametersPlots(Para,UserPara)
    Para2 = Dict()
    #= ------------------------------------------------------------
        Parameters for Plots
    ------------------------------------------------------------ =#
    Para2["ControlLabelsLatex"] = [latexstring("u_",kk) for kk=1:Para["nCon"]]
    Para2["ControlLabelsSimple"] = [string("u",kk) for kk=1:Para["nCon"]]
    Para2["ControlDistLabelsLatex"] = [latexstring("v_",kk) for kk=1:Para["nCon_dist"]]
    Para2["ControlDistLabelsSimple"] = [string("v",kk) for kk=1:Para["nCon_dist"]]
    Para2["StateLabelsLatex"] = [latexstring("x_",kk) for kk=1:Para["nStat"]]
    Para2["StateLabelsSimple"] = [string("x",kk) for kk=1:Para["nStat"]]
    Para2["StateDistLabelsLatex"] = [latexstring("y_",kk) for kk=1:Para["nStat_dist"]]
    Para2["StateDistLabelsSimple"] = [string("y",kk) for kk=1:Para["nStat_dist"]]
    Para2["CoStateLabelsLatex"] = [latexstring("\\lambda_",kk) for kk=1:Para["nStat"]]
    Para2["CoStateLabelsSimple"] = [string("lambda",kk) for kk=1:Para["nStat"]]
    Para2["CoStateDistLabelsLatex"] = [latexstring("\\xi_",kk) for kk=1:Para["nStat_dist"]]
    Para2["CoStateDistLabelsSimple"] = [string("xi_",kk) for kk=1:Para["nStat_dist"]]

    Para2["nVintagePlot"] = 10                      # Number of vintages plotted for realtime solutions
    Para2["PlotResultsStart"] = true
    Para2["PlotResultsIntermediate"] = true             # Indicator if intermediate results are plotted
    Para2["PlotResultsIntermediateFrequency"] = 40
    Para2["PlotResultsWaitForKey"] = false
    Para2["SavePlot"] = false                           # Indicator if plotted results should be exported and saved
    Para2["SavePlotPath"] = "Results/Plots/"            # Path for saved plots

    for kk in keys(Para2)
        if kk in keys(UserPara)
            Para[kk] = UserPara[kk]
        else
            Para[kk] = Para2[kk]
        end
    end
end


"""
    InitVariables(Para::Dict)

    Initialize control, state and costate variables for the Parameters specified
"""
function InitVariables(Para::Dict)
    Con = zeros(1,Para["nTime"],Para["nCon"])
    Con_dist = zeros(Para["nAge"],Para["nTime"],Para["nCon_dist"])
    Stat = zeros(1,Para["nTime"],Para["nStat"])
    Stat_dist = zeros(Para["nAge"],Para["nTime"],Para["nStat_dist"])
    Stat_agg = zeros(1,Para["nTime"],1)
    CoStat = zeros(1,Para["nTime"],Para["nStat"])
    CoStat_dist = zeros(Para["nAge"],Para["nTime"],Para["nStat_dist"])
    CoStat_agg = zeros(1,Para["nTime"],1)

    for kk in Para["Reduce1"]
        Con[:,:,kk] .= 0
    end
    for kk in Para["Reduce2"]
        Con_dist[:,:,kk] .= 0
    end

    return Con, Stat, Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,CoStat_agg
end


