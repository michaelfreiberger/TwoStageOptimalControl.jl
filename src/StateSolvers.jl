
"""
    state_PDE_solver(Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},Para::Dict)

    state_PDE_solver solves the dynamic system of state equations with given initial conditions and control variables.
"""
function state_PDE_solver(Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},Para::Dict)
    for ii=1:Para["nStat"]
        Stat[1,1,ii] = Para["InitStat"][ii]
    end

    #   Solve for the first stage states and consumption
    HeunODE(Con,Stat,Para)

    #   Initial values for second stage variables
    ShockTransition(Con,Stat,Stat_dist,Para)

    #   Loop in Time for vintage states jj
    if Para["ParallelComputing"] == true
        Threads.@threads for jj = 1:(Para["nTime"]-1)
            HeunPDE(jj,Con_dist,Stat_dist,Para)
        end
    else
        for jj = 1:(Para["nTime"]-1)
            HeunPDE(jj,Con_dist,Stat_dist,Para)
        end
    end

    #   Calculate aggregated variables
    Aggregate(Con_dist,Stat_dist,Stat_agg,Para)
end

"""
    HeunODE(Con::Array{Float64,3},Stat::Array{Float64,3},Para::Dict)

    HeunODEco solves the state dynamics in the first stage forward in time using the Heun-method for ordinary differential equations.
"""
function HeunODE(Con::Array{Float64,3},Stat::Array{Float64,3},Para::Dict)

    Dt = zeros(Para["nStat"])

    for ii = 1:(Para["nTime"]-1)
        #--------------------------------------------
        #   Heun method for ODE
        f_ODE(ii,Con,Stat,Dt,Para)
        Stat[1,ii+1,:] = Stat[1,ii,:] + Para["hstep"]*Dt
        
        f_ODE(ii+1,Con,Stat,Dt,Para)
        Stat[1,ii+1,:] = 0.5 * Stat[1,ii,:] + 0.5 * (Stat[1,ii+1,:] + Para["hstep"]*Dt)
        #   End Heun method for ODE
        #-----------------------------------------------
    end
end

"""
    RK4ODE(Con::Array{Float64,3},Stat::Array{Float64,3},Para::Dict)

    RK4ODE solves the state dynamics in the first stage forward in time using a 4th-order Runge-Kutta method for ordinary differential equations.
"""
function RK4ODE(Con::Array{Float64,3},Stat::Array{Float64,3},Para::Dict)
    c2 = 1/3
    c3 = 2/3
    c4 = 1.0

    a21 = 1/3
    a31 = -1/3
    a32 = 1.0
    a41 = 1.0
    a42 = -1.0
    a43 = 1.0

    b1 = 1/8
    b2 = 3/8
    b3 = 3/8
    b4 = 1/8

    for ii = 1:(Para["nTime"]-1)
        k1 = f_ODE_interstep(Para["tmesh"][ii],Con[1,ii,:],Stat,Para)
        k2 = f_ODE_interstep(Para["tmesh"][ii]+c2*Para["hstep"],
                        (1-c2)*Con[1,ii,:] + c2*Con[1,ii+1,:],
                        Stat[1,ii,:] + (a21*k1)*Para["hstep"],
                        Para)
        k3 = f_ODE_interstep(Para["tmesh"][ii]+c3*Para["hstep"],
                        (1-c3)*Con[1,ii,:] + c3*Con[1,ii+1,:],
                        Stat[1,ii,:] + (a31*k1 + a32*k2)*Para["hstep"],
                        Para)
        k4 = f_ODE_interstep(Para["tmesh"][ii]+c4*Para["hstep"],
                        (1-c4)*Con[1,ii,:] + c4*Con[1,ii+1,:],
                        Stat[1,ii,:] + (a41*k1 + a42*k2 + a43*k3)*Para["hstep"],
                        Para)

        Stat[1,ii+1,:] = Stat[1,ii,:] + Para["hstep"]*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
    end
end

"""
    HeunPDE(jj,Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Para::Dict)

    HeunODEco solves the state dynamics for vintage jj in the second stage forward in time using the Heun-method for ordinary differential equations.
"""
function HeunPDE(jj,Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Para::Dict)

    Dt = zeros(Para["nStat_dist"])

    for ii = jj:(Para["nTime"]-1)
        #--------------------------------------------
        #   Heun method for PDE
        f_PDE(ii,jj,Con_dist,Stat_dist,Dt,Para)
        Stat_dist[jj,ii+1,:] = Stat_dist[jj,ii,:] + Para["hstep"]*Dt
        
        f_PDE(ii+1,jj,Con_dist,Stat_dist,Dt,Para)
        Stat_dist[jj,ii+1,:] = 0.5 * Stat_dist[jj,ii,:] + 0.5 * (Stat_dist[jj,ii+1,:] + Para["hstep"]*Dt)
        #   End Heun method for PDE
        #-----------------------------------------------
    end
end

"""
    costate_PDE_solver(Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},
                            CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},CoStat_agg::Array{Float64,3},Para::Dict)

    costate_PDE_solver solves the dynamic system of co-state equations backwards with given endvalues for given controls and state variables.
    For the special case of endconstraints for the STATE-variables, we use the FOC to derive the respective value of CO-STATE variable.
"""
function costate_PDE_solver(Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},Stat_agg::Array{Float64,3},
                            CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},CoStat_agg::Array{Float64,3},Para::Dict)
    #Solves the PDE for the Costate equations for given controls u_M and u_H and state variables

    EndConstraintCoStat(Stat,CoStat,Para)
    EndConstraintCoStat_dist(Stat_dist,CoStat_dist,Para)

    if Para["ParallelComputing"] == true
        Threads.@threads for jj = 1:(Para["nTime"]-1)
            HeunPDEco(jj,Con_dist,Stat_dist,CoStat_dist,Para)
        end
    else
        for jj = 1:(Para["nTime"]-1)
            HeunPDEco(jj,Con_dist,Stat_dist,CoStat_dist,Para)
        end
    end
    HeunODEco(Con,Stat,CoStat,CoStat_dist,Para)
end

"""
    HeunODEco(Con::Array{Float64,3},Stat::Array{Float64,3},CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},Para::Dict)

    HeunODEco solves the costate dynamics in the first stage backwards in time using the Heun-method for ordinary differential equations.
"""
function HeunODEco(Con::Array{Float64,3},Stat::Array{Float64,3},CoStat::Array{Float64,3},CoStat_dist::Array{Float64,3},Para::Dict)
    Dt = zeros(Para["nStat"])

    #   Solves for non-distributed costates
    for ii = Para["nTime"]:(-1):2
        #--------------------------------------------
        #   Heun method for ODE for lambda
        f_ODE_co(ii,Con,Stat,CoStat,CoStat_dist,Dt,Para)
        CoStat[1,ii-1,:] =  CoStat[1,ii,:] - Para["hstep"]*Dt
        
        f_ODE_co(ii-1,Con,Stat,CoStat,CoStat_dist,Dt,Para)
        CoStat[1,ii-1,:] = 0.5*CoStat[1,ii,:] + 0.5*(CoStat[1,ii-1,:] - Para["hstep"]*Dt)
        #   End Heun method for ODE for lambda
        #-----------------------------------------------
    end
end

"""
    HeunPDEco(jj,Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},CoStat_dist::Array{Float64,3},Para::Dict)

    HeunPDEco solves the costate dynamic for vintage jj in the second stage backwards in time using the Heun-method for ordinary differential equations.
"""
function HeunPDEco(jj,Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3},CoStat_dist::Array{Float64,3},Para::Dict)
    
    Dt = zeros(Para["nStat_dist"])
    #   Solve costates for vintage jj
    for ii = Para["nTime"]:(-1):(jj+1)
        
        #---------------------------------------
        #   Heun method for xi
        f_PDE_co(ii,jj,Con_dist,Stat_dist,CoStat_dist,Dt,Para)
        CoStat_dist[jj,ii-1,:] = CoStat_dist[jj,ii,:] - Para["hstep"]*Dt
        
        f_PDE_co(ii-1,jj,Con_dist,Stat_dist,CoStat_dist,Dt,Para)
        CoStat_dist[jj,ii-1,:] = 0.5*CoStat_dist[jj,ii,:] + 0.5*(CoStat_dist[jj,ii-1,:] - Para["hstep"]*Dt)
        #   End Heun method for xi
        #----------------------------------------
    end
end

