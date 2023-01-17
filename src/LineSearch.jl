
"""
    LineSearch(Con::Array{Float64,3}, Stat::Array{Float64,3}, Con_dist::Array{Float64,3},  Stat_dist::Array{Float64,3}, Stat_agg::Array{Float64,3}, 
                    Dir_u::Array{Float64,3}, Dir_u_dist::Array{Float64,3}, ObjValue::Float64, Step, Para::Dict,
                    Con_new::Array{Float64,3}, Stat_new::Array{Float64,3}, Con_dist_new::Array{Float64,3}, Stat_dist_new::Array{Float64,3}, Stat_agg_new::Array{Float64,3},
                    Con_best::Array{Float64,3}, Stat_best::Array{Float64,3}, Con_dist_best::Array{Float64,3}, Stat_dist_best::Array{Float64,3}, Stat_agg_best::Array{Float64,3})

    Searches for an improvement along the gradient using a quadratic approximation for the objective function along the gradient.
"""
function LineSearch(Con::Array{Float64,3}, Stat::Array{Float64,3}, Con_dist::Array{Float64,3},  Stat_dist::Array{Float64,3}, Stat_agg::Array{Float64,3}, 
                    Dir_u::Array{Float64,3}, Dir_u_dist::Array{Float64,3}, ObjValue::Float64, Step, Para::Dict,
                    Con_new::Array{Float64,3}, Stat_new::Array{Float64,3}, Con_dist_new::Array{Float64,3}, Stat_dist_new::Array{Float64,3}, Stat_agg_new::Array{Float64,3},
                    Con_best::Array{Float64,3}, Stat_best::Array{Float64,3}, Con_dist_best::Array{Float64,3}, Stat_dist_best::Array{Float64,3}, Stat_agg_best::Array{Float64,3})
            # Perform line search for new optimal control

    impr = 0
    Para["LineIter"] = 0

    ObjValue_new = 0
    ObjValue_best = 0

    while (impr==0) && (Para["LineIter"]<=Para["MaxLineIter"]) && (Step > Para["stepLowBound"])

        Step = min(Step,Para["UpperLineStep"])

        ObjValue_new, Step = Update(Con, Con_dist, Dir_u, Dir_u_dist,Step,Para,
                            Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new)

        y0 = ObjValue
        y1 = ObjValue_new
        x1 = Step

        if y1 > y0
            impr = 1

            # Save currently best values
            ObjValue_best, step_best = AssignBest(ObjValue_new, Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new, Step, Con_best, Stat_best, Con_dist_best, Stat_dist_best, Stat_agg_best)

            step_increase = Para["LineSearchStepIncrease"]
            Step = (1+step_increase)*Step
            x2 = Step

            ObjValue_new, Step = Update(Con, Con_dist, Dir_u, Dir_u_dist,Step,Para,
                                    Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new)
            y2 = ObjValue_new

            if y2 > y1
                # Save currently best values
                ObjValue_best, step_best = AssignBest(ObjValue_new, Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new, Step, Con_best, Stat_best, Con_dist_best, Stat_dist_best, Stat_agg_best)

                if x1*y2 + x2*y0 - x1*y0 - x2*y1 <= 0       # Check if curve is convex
                    # Curve is convex --> second solution is the best
                else
                    # Curve is concave --> stepsize for maximum of quadratic approximation
                    st = ( (x2^2-x1^2)*y0 + y2*x1^2 - y1*x2^2 ) / ( 2*(y2*x1 - y1*x2 + y0*(x2-x1)) )
                    st = min(st,Para["UpperLineStep"])

                    ObjValue_new, st = Update(Con, Con_dist, Dir_u, Dir_u_dist,Step,Para,
                                                Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new)

                    if ObjValue_new > y2
                        # Save currently best values
                        ObjValue_best, step_best = AssignBest(ObjValue_new, Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new, st, Con_best, Stat_best, Con_dist_best, Stat_dist_best, Stat_agg_best)
                    end
                end
            else    #(y2 <= y1)
                # Curve is concave --> stepsize for maximum of quadratic approximation
                st = ( (x2^2-x1^2)*y0 + y2*x1^2 - y1*x2^2 ) / ( 2*(y2*x1 - y1*x2 + y0*(x2-x1)) )
                st = min(st,Para["UpperLineStep"])

                ObjValue_new, st = Update(Con, Con_dist, Dir_u, Dir_u_dist,Step,Para,
                                        Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new)

                if ObjValue_new > y1
                    # Save currently best values
                    ObjValue_best, step_best = AssignBest(ObjValue_new, Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new, st, Con_best, Stat_best, Con_dist_best, Stat_dist_best, Stat_agg_best)
                end
            end
        else     # (y1 <= y0)
            Step = 0.5*Step
            x2 = Step

            ObjValue_new, Step = Update(Con, Con_dist, Dir_u, Dir_u_dist,Step,Para,
                                    Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new)
            y2 = ObjValue_new

            if y2 > y0
                impr = 1
                # Save currently best values
                ObjValue_best, step_best = AssignBest(ObjValue_new, Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new, Step, Con_best, Stat_best, Con_dist_best, Stat_dist_best, Stat_agg_best)

                # Curve is concave --> stepsize for maximum of quadratic approximation
                st = ( (x2^2-x1^2)*y0 + y2*x1^2 - y1*x2^2 ) / ( 2*(y2*x1 - y1*x2 + y0*(x2-x1)) )
                st = min(st,Para["UpperLineStep"])

                # New controls + new states + new objective function()
                ObjValue_new, st = Update(Con, Con_dist, Dir_u, Dir_u_dist,Step,Para,
                                        Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new)

                if ObjValue_new > y2
                    # Save currently best values
                    ObjValue_best, step_best = AssignBest(ObjValue_new, Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new, st, Con_best, Stat_best, Con_dist_best, Stat_dist_best, Stat_agg_best)
                end
            elseif y2 > 0.5*(y0+y1)     # && y2 < (3*y0+y1)/4
                st = ( (x2^2-x1^2)*y0 + y2*x1^2-y1*x2^2 ) / ( 2*(y2*x1 - y1*x2 + y0*(x2-x1)) )
                st = min(st,Para["UpperLineStep"])
                if st > 0
                    ObjValue_new, st = Update(Con, Con_dist, Dir_u, Dir_u_dist,Step,Para,
                                            Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new)

                    if ObjValue_new >= y0
                        impr = 1
                        # Save currently best values
                        ObjValue_best, step_best = AssignBest(ObjValue_new, Con_new, Stat_new, Con_dist_new, Stat_dist_new, Stat_agg_new, st, Con_best, Stat_best, Con_dist_best, Stat_dist_best, Stat_agg_best)
                    end
                end
            end
        end
        if Para["SearchDISP"] == true
            if impr==1
                printSuccess(Para, ObjValue_best, step_best)
            else
                printNoSuccess(Para, ObjValue, Step)
            end
        end
    end

    if impr == 0
    # If there is no improvement, return original values
    step_best = Step
    # only step-size gets adjusted --> new try from the beginning with the new step-size
    else
    # Assign the best values to the basic variables without _best or _new extension
    ObjValue, Step = AssignBest(ObjValue_best, Con_best, Stat_best, Con_dist_best, Stat_dist_best, Stat_agg_best, step_best, Con, Stat, Con_dist, Stat_dist, Stat_agg)
    end

    return ObjValue,step_best,impr
end


"""
    Update(Con::Array{Float64,3},Con_dist::Array{Float64,3},Dir_u::Array{Float64,3},Dir_u_dist::Array{Float64,3},Step,Para::Dict,
                Con_new::Array{Float64,3},Stat_new::Array{Float64,3},Con_dist_new::Array{Float64,3},Stat_dist_new::Array{Float64,3}, Stat_agg_new::Array{Float64,3})

    Update calculates the new values for Controls `Con_new`, States `Stat_new` and Objectiv value given start values for the Controls `Con`, the search direction `Dir_u` and the step-size in the direction `step`
"""
function Update(Con::Array{Float64,3},Con_dist::Array{Float64,3},Dir_u::Array{Float64,3},Dir_u_dist::Array{Float64,3},Step,Para::Dict,
                Con_new::Array{Float64,3},Stat_new::Array{Float64,3},Con_dist_new::Array{Float64,3},Stat_dist_new::Array{Float64,3}, Stat_agg_new::Array{Float64,3})

    Step = min(Step,Para["UpperLineStep"])

    Con_new .= Con + Step*Dir_u
    Con_dist_new .= Con_dist + Step*Dir_u_dist
    ConMapping(Con_new,Para)
    for ii = 1:Para["nAge"]
    ConMapping_dist(Con_dist_new,ii,Para)
    end

    state_PDE_solver(Con_new,Stat_new,Con_dist_new, Stat_dist_new, Stat_agg_new, Para)

    ObjValue_new = ObjectValue(Con_new,Stat_new,Con_dist_new,Stat_dist_new,Stat_agg_new,Para)
    Para["LineIter"] = Para["LineIter"] + 1

    return ObjValue_new, Step
end

"""
    AssignBest(ObjValue::Float64,Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3}, Stat_agg::Array{Float64,3}, 
                    step, 
                    Con_best::Array{Float64,3}, Stat_best::Array{Float64,3}, Con_dist_best::Array{Float64,3},  Stat_dist_best::Array{Float64,3},  Stat_agg_best::Array{Float64,3})

    Assign the current variables to the "best" variables
"""
function AssignBest(ObjValue::Float64,Con::Array{Float64,3},Stat::Array{Float64,3},Con_dist::Array{Float64,3},Stat_dist::Array{Float64,3}, Stat_agg::Array{Float64,3}, 
                    step, 
                    Con_best::Array{Float64,3}, Stat_best::Array{Float64,3}, Con_dist_best::Array{Float64,3},  Stat_dist_best::Array{Float64,3},  Stat_agg_best::Array{Float64,3})

    Con_best .= Con
    Stat_best .= Stat
    Con_dist_best .= Con_dist
    Stat_dist_best .= Stat_dist
    Stat_agg_best .= Stat_agg

    return ObjValue, step
end
