

"""
    ObjectValue( Con::Array{Float64,3}, Stat::Array{Float64,3}, Con_dist::Array{Float64,3}, Stat_dist::Array{Float64,3},  Stat_agg::Array{Float64,3}, Para::Dict )

    Calculate the aggregated objective value of the optimal control problem for the current solutions of control and state variables
"""
function ObjectValue( Con::Array{Float64,3}, Stat::Array{Float64,3}, Con_dist::Array{Float64,3}, Stat_dist::Array{Float64,3},  Stat_agg::Array{Float64,3}, Para::Dict )
    F = zeros(Para["nTime"])
    for jj = 1:Para["nTime"]
        F[jj] = exp(-Para["rho"]*Para["tmesh"][jj]) * 
                (ObjectiveIntegrand(Con[1,jj,:],Stat[1,jj,:],Para["tmesh"][jj],Para) + Stat_agg[1,jj,1])
    end
    return Integral(Para["nTime"],Para["hstep"],F)
end

"""
    Aggregate(Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},Para::Dict )

    Calculate the aggregated variables over all vintage in the second stage
"""
function Aggregate(Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},Para::Dict )
    for ii = 1:Para["nTime"]
        F = zeros(ii)
        for jj = 1:ii
            F[jj] = AggregationFunction(Con_dist[jj,ii,:],Stat_dist[jj,ii,:],Para["tmesh"][ii],Para["tmesh"][jj],Para)
        end
        Stat_agg[1,ii,1] = Integral(ii,Para["hstep"],F)
    end
end

"""
    ShockTransition(Con::Array{Float64,3},Stat::Array{Float64,3},Stat_dist::Array{Float64,3},Para::Dict )

    Calculate the starting values for all vintages given the state profiles of the first stage.
"""
function ShockTransition(Con::Array{Float64,3},Stat::Array{Float64,3},Stat_dist::Array{Float64,3},Para::Dict )
    for ii = 1:Para["nTime"]
        Stat_dist[ii,ii,:] = Shock(Con[1,ii,:],Stat[1,ii,:],Para["tmesh"][ii],Para)
    end
end

"""
    f_ODE( tt::Int64, Con::Array{Float64,3}, Stat::Array{Float64,3}, Dt::Array{Float64,1}, Para::Dict )

    Right hand side of the ODE for the state variables in the first stage at time index tt.
"""
function f_ODE( tt::Int64, Con::Array{Float64,3}, Stat::Array{Float64,3}, Dt::Array{Float64,1}, Para::Dict )
    Dt .= StateDynamic_1(Con[1,tt,:],Stat[1,tt,:],Para["tmesh"][tt],Para)
end

"""
    f_ODE_interstep( t::Float64, Con, Stat, Para::Dict)

    Right hand side of the ODE for the state variables in the first stage. Only used for the 4th order RK-method where the Control variables have to be interpolated.
"""
function f_ODE_interstep( t::Float64, Con, Stat, Para::Dict)
    k = StateDynamic_1(Con,Stat,t,Para)
    return k
end

"""
    f_PDE( tt::Int64,ss::Int64, Con_dist::Array{Float64,3}, Stat_dist::Array{Float64,3}, Dt::Array{Float64,1},Para::Dict )

    Right hand side of the PDE for the state variables in the second stage at time index tt and vintage ss.
"""
function f_PDE( tt::Int64,ss::Int64, Con_dist::Array{Float64,3}, Stat_dist::Array{Float64,3}, Dt::Array{Float64,1},Para::Dict )
    Dt .= StateDynamic_2(Con_dist[ss,tt,:],Stat_dist[ss,tt,:],Para["tmesh"][tt],Para["tmesh"][ss],Para)
end

"""
    f_ODE_co( tt::Int64,Con::Array{Float64,3},Stat::Array{Float64,3},CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},Dt::Array{Float64,1},Para::Dict )

    Right hand side of the PDE for the costate variables in the first stage at time index tt.
"""
function f_ODE_co( tt::Int64,Con::Array{Float64,3},Stat::Array{Float64,3},CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},Dt::Array{Float64,1},Para::Dict )
    Dt .=  Para["rho"]*CoStat[1,tt,:] - ForwardDiff.gradient(X->Hamiltonian(Con[1,tt,:],X,CoStat[1,tt,:],CoStat_dist[tt,tt,:],Para["tmesh"][tt],Para),Stat[1,tt,:])
end

"""
    f_PDE_co( tt::Int64,ss::Int64,Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},CoStat_dist::Array{Float64,3},Dt::Array{Float64,1},Para::Dict )

    Right hand side of the PDE for the costate variables in the second stage at time index tt and vintage ss.
"""
function f_PDE_co( tt::Int64,ss::Int64,Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},CoStat_dist::Array{Float64,3},Dt::Array{Float64,1},Para::Dict )
    Dt .= Para["rho"]*CoStat_dist[ss,tt,:] - ForwardDiff.gradient(X->Hamiltonian_dist(Con_dist[ss,tt,:],X,CoStat_dist[ss,tt,:],Para["tmesh"][tt],Para["tmesh"][ss],Para),Stat_dist[ss,tt,:])
end

"""
    Hamiltonian(Con, Stat, CoStat, CoStat_dist, t::Float64, Para::Dict )

    Definition of the Hamiltonian with all terms relevant for the calculations including first stage variables.
"""
function Hamiltonian(Con, Stat, CoStat, CoStat_dist, t::Float64, Para::Dict )
    return ObjectiveIntegrand(Con,Stat,t,Para) + 
            + dot(CoStat,StateDynamic_1(Con,Stat,t,Para)) + 
            + dot(CoStat_dist,Shock(Con,Stat,t,Para))
end

"""
    Hamiltonian_dist(Con_dist, Stat_dist, CoStat_dist, t::Float64, s::Float64, Para::Dict )

    Definition of the Hamiltonian with all terms relevant for the calculations including second stage variables.
"""
function Hamiltonian_dist(Con_dist, Stat_dist, CoStat_dist, t::Float64, s::Float64, Para::Dict )
    return dot(CoStat_dist,StateDynamic_2(Con_dist,Stat_dist,t,s,Para)) +
            + AggregationFunction(Con_dist,Stat_dist, t,s, Para)
end

"""
    GradHamiltonian(Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},dHam::Array{Float64,3},dHam_dist::Array{Float64,3},Para::Dict)

    Calculation of the gradient of the Hamiltonian for both stages.
"""
function GradHamiltonian(Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},
    CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},dHam::Array{Float64,3},dHam_dist::Array{Float64,3},Para::Dict)
    
    dHam_dist .= 0
    dHam .= 0

    if Para["ParallelComputing"] == true
        Threads.@threads for jj = 1:Para["nTime"]    
            # Calculate Gradient for first stage
            dHam[1,jj,:] = ForwardDiff.gradient(X->Hamiltonian(X,Stat[1,jj,:],CoStat[1,jj,:],CoStat_dist[jj,jj,:],Para["tmesh"][jj],Para),Con[1,jj,:])
            # Calculate Gradient for vintage jj
            for ii = jj:Para["nTime"]
                dHam_dist[jj,ii,:] = ForwardDiff.gradient(X->Hamiltonian_dist(X,Stat_dist[jj,ii,:],CoStat_dist[jj,ii,:],Para["tmesh"][ii],Para["tmesh"][jj],Para),Con_dist[jj,ii,:])
            end
        end
    else
        for jj = 1:Para["nTime"]    
            # Calculate Gradient for first stage
            dHam[1,jj,:] = ForwardDiff.gradient(X->Hamiltonian(X,Stat[1,jj,:],CoStat[1,jj,:],CoStat_dist[jj,jj,:],Para["tmesh"][jj],Para),Con[1,jj,:])
            # Calculate Gradient for vintage jj
            for ii = jj:Para["nTime"]
                dHam_dist[jj,ii,:] = ForwardDiff.gradient(X->Hamiltonian_dist(X,Stat_dist[jj,ii,:],CoStat_dist[jj,ii,:],Para["tmesh"][ii],Para["tmesh"][jj],Para),Con_dist[jj,ii,:])
            end
        end
    end
end

"""
    NewDirection(Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},dHam::Array{Float64,3},dHam_dist::Array{Float64,3},Para::Dict)

    Adjustment of the gradient of the Hamiltonian depending on the optimisation type chosen. Available types are:
        - "Gradient"
            Uses the gradient without adjustments
        - "ProbAdjust"
            Adjusts the gradient for the weighting of the vintage in the second stage and the survival probability in the first stage
        - "Newton-Raphson"
            Adjusts the gradient as described in the Newton-Raphson method using the Hessian of the Hamiltonian
"""
function NewDirection(Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},dHam::Array{Float64,3},dHam_dist::Array{Float64,3},Para::Dict)
    eps = 1e-6
    if Para["OptiType"] == "Gradient"
        dHam .= dHam
        dHam_dist .= dHam_dist
    elseif Para["OptiType"] == "ProbAdjust"
        for jj = 1:Para["nTime"]
            for ii = jj:Para["nTime"]
                dHam_dist[jj,ii,:] .= dHam_dist[jj,ii,:] ./ max(eps,Stat_dist[jj,ii,Para["ProbIndex"]])
            end
        end
        for ii = 1:Para["nTime"]
            dHam[1,ii,:] .= dHam[1,ii,:] ./ max(eps,Stat[1,ii,Para["ProbIndex"]])
        end
    elseif Para["OptiType"] == "Newton-Raphson"
        Hessian(ii,jj) = ForwardDiff.hessian(X->Hamiltonian_dist(X,Stat_dist[jj,ii,:],CoStat_dist[jj,ii,:],Para["tmesh"][ii],Para["tmesh"][jj],Para),Con_dist[jj,ii,:])
        for jj = 1:Para["nTime"]
            for ii = jj:Para["nTime"]
                dHam_dist[jj,ii,:] .= - inv(Hessian(ii,jj))*dHam_dist[jj,ii,:] 
            end
        end
        Hessian(ii) = ForwardDiff.hessian(X->Hamiltonian(X,Stat[1,ii,:],CoStat[1,ii,:],CoStat_dist[ii,ii,:],Para["tmesh"][ii],Para),Con[1,ii,:])
        for ii = 1:Para["nTime"]
            dHam[1,ii,:] .= -inv(Hessian(ii))*dHam[1,ii,:]
        end
    else
        Para["OptiType"] = "Gradient"
    end

    #=----------------------------
        Smoothing
    ----------------------------=#
    if Para["GradSmooth"] == true
        GradSmooth(dHam_dist,dHam,Para)
    end

    #=-------------------------------------------------
        Reducing the model
    -------------------------------------------------=#
    for kk in Para["Reduce1"]
        dHam[:,:,kk] .= 0
    end
    for kk in Para["Reduce2"]
        dHam_dist[:,:,kk] .= 0
    end
end

"""
    ConMapping( Con::Array{Float64,3}, Para )

    Map control variables in the first stage into the feasible region described by the lower and upper bounds of the controls.
"""
function ConMapping( Con::Array{Float64,3}, Para )
    for kk = 1:Para["nCon"]
        if !(kk in Para["Reduce1"])
            Con[1,:,kk] .= max.(Con[1,:,kk],Para["ConMin"][kk])
            Con[1,:,kk] .= min.(Con[1,:,kk],Para["ConMax"][kk])
        end
    end
end

"""
    ConMapping_dist( Con_dist::Array{Float64,3}, ii, Para )

    Map control variables in the second stage into the feasible region described by the lower and upper bounds of the controls.
"""
function ConMapping_dist( Con_dist::Array{Float64,3}, ii, Para )
    for kk = 1:Para["nCon_dist"]
        if !(kk in Para["Reduce2"])
            Con_dist[ii,ii:end,kk] = max.(Con_dist[ii,ii:end,kk],Para["Con_distMin"][kk])
            Con_dist[ii,ii:end,kk] = min.(Con_dist[ii,ii:end,kk],Para["Con_distMax"][kk])
        end
    end
end

"""
    EndConstraintCoStat(Stat::Array{Float64,3},CoStat::Array{Float64,3},Para::Dict)

    Define the endconstraint for the first stage costate variables using the derivative of the salvage-function.
"""
function EndConstraintCoStat(Stat::Array{Float64,3},CoStat::Array{Float64,3},Para::Dict)
    CoStat[1,Para["nTime"],:] = ForwardDiff.gradient(X->SalvageFunction_1(X,Para),Stat[1,Para["nTime"],:])
end

"""
    EndConstraintCoStat_dist(Stat_dist::Array{Float64,3},CoStat_dist::Array{Float64,3},Para::Dict)

    Define the endconstraint for the second stage costate variables using the derivative of the salvage-function.
"""
function EndConstraintCoStat_dist(Stat_dist::Array{Float64,3},CoStat_dist::Array{Float64,3},Para::Dict)
    for ii = 1:Para["nAge"]
        CoStat_dist[ii,Para["nTime"],:] = ForwardDiff.gradient(X->SalvageFunction_2(X,Para),Stat_dist[ii,Para["nTime"],:])
    end
end


