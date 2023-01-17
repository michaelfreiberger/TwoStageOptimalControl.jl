

"""
    TwoStageOptimisation(;Results=Dict(),UserParameters=Dict(),
                    ObjectiveIntegrand2 = ObjectiveIntegrand,AggregationFunction2 = AggregationFunction,
                    StateDynamic_1_2 = StateDynamic_1, StateDynamic_2_2 = StateDynamic_2, Shock2 = Shock,
                    SalvageFunction_1_2 = SalvageFunction_1, SalvageFunction_2_2 = SalvageFunction_2)

    
"""
function TwoStageOptimisation(;Results=Dict(),UserParameters=Dict(),
                    ObjectiveIntegrand2 = ObjectiveIntegrand,AggregationFunction2 = AggregationFunction,
                    StateDynamic_1_2 = StateDynamic_1, StateDynamic_2_2 = StateDynamic_2, Shock2 = Shock,
                    SalvageFunction_1_2 = SalvageFunction_1, SalvageFunction_2_2 = SalvageFunction_2)
    
    
    # Set up dictionary with all parameters either from the base values or user-specific values
    Para = Dict()
    ParametersBasics(Para,UserParameters)
    ParametersGrids(Para,UserParameters)
    ParametersControls(Para,UserParameters)
    ParametersAlgorithm(Para,UserParameters)
    ParametersInitStates(Para,UserParameters)
    ParametersPlots(Para,UserParameters)
    ParametersFunctions(Para,UserParameters)
    for kk in keys(UserParameters)
        if !(kk in keys(Para))
            Para[kk] = UserParameters[kk]
        end
    end

    # Initialize the control, state and costate variables
    Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,CoStat_agg = InitVariables(Para)

    #-----------------------------------------------------------------------------------------------
    # Check wether all the function supplied by the user have the right dimensions
    SystemCheck = 1

    if !(ObjectiveIntegrand2(Con[1,1,:],Stat[1,1,:],0.0,Para) isa Number)
        global ObjectiveIntegrand = (Con,Stat, t::Float64, Para::Dict) -> 0.0
        println("Warning: Objective Function is not a number! Objective Function is set to 0")
        SystemCheck = 0
    else
        global ObjectiveIntegrand = ObjectiveIntegrand2
    end
    if !(AggregationFunction2(Con[1,1,:],Stat[1,1,:],0.0,0.0,Para) isa Number)
        global AggregationFunction = (Con,Stat, t::Float64,s::Float64, Para::Dict) -> 0.0
        println("Warning: Aggregation Function is not a number! Aggregation Function is set to 0")
        SystemCheck = 0
    else
        global AggregationFunction = AggregationFunction2
    end
    if size(StateDynamic_1_2(Con[1,1,:],Stat[1,1,:],0.0,Para),1) != Para["nStat"]
        global StateDynamic_1 = (Con,Stat, t::Float64, Para::Dict) -> zeros(Para["nStat"])
        println("Warning: Dimensions of State Dynamics (Stage 1) do not match! Dynamics are set to 0")
        SystemCheck = 0
    else
        global StateDynamic_1 = StateDynamic_1_2
    end
    if size(StateDynamic_2_2(Con_dist[1,1,:],Stat_dist[1,1,:],0.0,0.0,Para),1) != Para["nStat_dist"]
        global StateDynamic_2 = (Con_dist,Stat_dist, t::Float64, s::Float64, Para::Dict) -> zeros(Para["nStat_dist"])
        println("Warning: Dimensions of State Dynamics (Stage 2) do not match! Dynamics are set to 0")
        SystemCheck = 0
    else
        global StateDynamic_2 = StateDynamic_2_2
    end
    if size(Shock2(Con[1,1,:],Stat[1,1,:],0.0,Para),1) != Para["nStat_dist"]
        global Shock = (Con, Stat, t::Float64, Para::Dict) -> zeros(Para["nStat_dist"])
        println("Warning: Shock Transition dimensions do not math! Transition is set to 0!")
        SystemCheck = 0
    else
        global Shock = Shock2
    end
    if !(SalvageFunction_1_2(Stat[1,1,:],Para) isa Number)
        global SalvageFunction_1 = (Stat,Para::Dict) -> 0.0
        println("Warning: Salvage Function (Stage 1) is not a number! Salvage Function is set to 0")
        SystemCheck = 0
    else
        global SalvageFunction_1 = SalvageFunction_1_2
    end
    if !(SalvageFunction_2_2(Stat_dist[1,1,:],Dict()) isa Number)
        global SalvageFunction_2 = (Stat,Para::Dict) -> 0.0
        println("Warning: Salvage Function (Stage 2) is not a number! Salvage Function is set to 0")
        SystemCheck = 0
    else
        global SalvageFunction_2 = SalvageFunction_2_2
    end
    
    # Return error message if some dimensions do not match
    if SystemCheck == 0
        println("Calculations are terminated due to dimension mismatches!")
        return
    end
   

    #-----------------------------------------------------------------------------------------------
    # Define the order of the integration method for the objective value and 
    # the aggregated variables
    if Para["IntegrationOrder"] == 1
        global Integral = integ1
    elseif Para["IntegrationOrder"] == 2
        global Integral = integ2
    elseif Para["IntegrationOrder"] == 4
        global Integral = integ4
    else
        global Integral = integ4
    end


    #-----------------------------------------------------------------------------------------------
    # Load (interpolate) the supplied controls in the results dictionary
    if Para["LoadInits"] == true
        if all(in.(["Con", "Con_dist"],(keys(Results),)))
            Con, Con_dist = LoadVariables(Para,Results)
        else
            println("Warning: Initial Values cannot be loaded. Con and/or Con_dist are not elemet of the dictionary.")
        end
    end



    #------------------------------------------------------------------------
    #   Start Optimisation loops
    Para["HstepIter"] = 1
    Step = Para["InitLineStep"]

    # Stop calculations if the current time-step size is smaller or equal to the lower bound specified
    if Para["hstep"] > Para["hLowBound"]
        StopOuterIteration = 0
    else
        StopOuterIteration = 1
    end
    while StopOuterIteration == 0

        # Adjust parameters to the new time-step size
        ParaAdjust(Para["hstep"],Para)
        
        # Interpolate controls after hstep reduction
        Con, Con_dist = ConInterpol(Con,Con_dist,Para)

        # Initialize the gradients and other auxiliary variables
        dHam,Stat,dHam_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,CoStat_agg = InitVariables(Para)
        Con_new,Stat_new,Con_dist_new,Stat_dist_new,Stat_agg_new = InitVariables(Para)
        Con_best,Stat_best,Con_dist_best,Stat_dist_best,Stat_agg_best = InitVariables(Para)

        # Calculate the state variables related to the current control profiles
        state_PDE_solver(Con, Stat, Con_dist, Stat_dist, Stat_agg, Para)
        # Calculate the corresponding objective value
        ObjValue = ObjectValue(Con,Stat,Con_dist,Stat_dist,Stat_agg,Para)
        # Calculate the profiles of the costate variables
        costate_PDE_solver(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,CoStat_agg,Para)
        # Calculate the gradient of the Hamiltonian for the current control values
        GradHamiltonian(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,dHam,dHam_dist,Para)
        # Adjust gradient according to user specifications
        NewDirection(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,dHam,dHam_dist,Para)
        #Assign the current profiles to the Results dictionary
        AssignResults(Results, Con, Stat, Con_dist, Stat_dist, Stat_agg, CoStat, CoStat_dist, dHam, dHam_dist, Para)
        
        if Para["PlotResultsStart"]
            # Plot all profiles
            PlotResults(Results)
            # Wait for user key if specified
            if Para["PlotResultsWaitForKey"]
                display("Wait for key")
                readline()
            end
        end


        impr = 1
        while Para["OptiIter"] <= Para["MaxOptiIter"] && Step > Para["stepLowBound"] && impr==1  # basic iteration loop

            ObjValue = ObjectValue(Con,Stat,Con_dist,Stat_dist,Stat_agg,Para)
            costate_PDE_solver(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,CoStat_agg,Para)
            GradHamiltonian(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,dHam,dHam_dist,Para)
            NewDirection(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,dHam,dHam_dist,Para)

            # Linesearch along the gradient for a better set of controls
            ObjValue, Step, impr = LineSearch(Con,Stat,Con_dist,Stat_dist,Stat_agg,dHam,dHam_dist,ObjValue,Step,Para,
                                              Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new,
                                              Con_best, Stat_best, Con_dist_best, Stat_dist_best, Stat_agg_best)

            Para["OptiIter"] = Para["OptiIter"] + 1
            
            # Plot the current results every Para["PlotResultsIntermediateFrequency"] iterations
            if Para["PlotResultsIntermediate"] == true && Para["OptiIter"]%Para["PlotResultsIntermediateFrequency"] == 0
                GradHamiltonian(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,dHam,dHam_dist,Para)
                NewDirection(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,dHam,dHam_dist,Para)
                AssignResults(Results, Con, Stat, Con_dist, Stat_dist, Stat_agg, CoStat, CoStat_dist, dHam, dHam_dist, Para)
                PlotResults(Results)
                if Para["PlotResultsWaitForKey"]
                    display("Wait for key")
                    readline()
                end
            end
            
            # Smooth the profiles of the control variables if specified by the user
            if Para["OptiIter"]%50 == 0 && Para["ConSmooth"]
                AssignResults(Results, Con_dist, Con, Stat_dist, Stat, Stat_agg, CoStat_dist, CoStat, dHam_dist, dHam, Para)
                for ii = 1:10
                    ConMapping(Con_new,Para)
                    for ii = 1:Para["nAge"]
                        ConMapping_dist(Con_dist_new,ii,Para)
                    end
                    ConSmooth(Con_dist_new,Con_new,Para)
                end
                AssignResults(Results, Con_dist, Con, Stat_dist, Stat, Stat_agg, CoStat_dist, CoStat, dHam_dist, dHam, Para)
            end
        end

        AssignResults(Results, Con, Stat, Con_dist, Stat_dist, Stat_agg, CoStat, CoStat_dist, dHam, dHam_dist, Para)
        if Para["PlotResultsIntermediate"]
            PlotResults(Results)
        end

        # Decrease/Half the time-step size and reset the iteration counters
        if Para["hstep"]/2 >= Para["hLowBound"]
            Para["hstep"] = Para["hstep"] * 0.5
            Para["HstepIter"] = Para["HstepIter"] + 1
            Para["MaxOptiIter"] = floor(Int,Para["MaxOptiIter"] * 2)
            Para["UpperLineStep"] = Para["UpperLineStep"] * 1
            Step = max(Step,Para["InitLineStep"])
            StopOuterIteration = 0
            # if the smallest step-size is reached, use the original gradient for the linesearch without adjustment
            if Para["hstep"]/2 <= Para["hLowBound"] 
                Para["OptiType"] = "Gradient"
            end
        else
            StopOuterIteration = 1
        end
    end
    
    #----------------------------------------------------------------------------------------------------------
    # Calculate all final variables
    dHam,Stat,dHam_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,CoStat_agg = InitVariables(Para)
    state_PDE_solver(Con, Stat, Con_dist, Stat_dist, Stat_agg, Para)
    ObjValue = ObjectValue(Con,Stat,Con_dist,Stat_dist,Stat_agg,Para)
    costate_PDE_solver(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,CoStat_agg,Para)
    GradHamiltonian(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,dHam,dHam_dist,Para)
    Para["OptiType"] = "Gradient"
    NewDirection(Con,Stat,Con_dist,Stat_dist,Stat_agg,CoStat,CoStat_dist,dHam,dHam_dist,Para)

    AssignResults(Results, Con, Stat, Con_dist, Stat_dist, Stat_agg, CoStat, CoStat_dist, dHam, dHam_dist, Para)
    PlotResults(Results;SavePlot = Para["SavePlot"])
    return Results
end
