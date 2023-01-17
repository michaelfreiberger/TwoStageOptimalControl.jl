"""
    integ4(n::Int64,h::Float64,F)

    integ4 performs 4th order integration (Simpson-method) method for datavector `F` with `n` points with a distance of `h` between each of them.
    If n > 1 and F = 1, F will be interpreted as a constant function over the n points.
"""
function integ4(n::Int64,h::Float64,F)
   if n != length(F) && length(F) == 1
      # F is assumed to be constant
      F = F*ones(n)
   elseif n > length(F)
       println("n does not math number of data points (integ4)")
       println("n = ",n)
       println("F = ",F)
       return NaN
   end
   integ = 0
   if n < 1
       return 0
   elseif n==1
       integ = 0
   elseif n==2
       # trapezoidal rule
       integ = h*0.5*(F[1] + F[2])
   elseif n==3
       integ = 1/3 * (F[1] + 4*F[2] + F[3])
   elseif n==4
       integ = 3/8 * (F[1] + 3*F[2] + 3*F[3] + F[4])
   else
       for ii = 3:2:(n-2)
            integ = integ + 1/3*(F[ii-2] + 4*F[ii-1] + F[ii])
       end
       if isodd(n)
            integ = integ + 1/3*(F[n-2] + 4*F[n-1] + F[n])
       elseif iseven(n)
            integ = integ + 3/8 * (F[n-3] + 3*F[n-2] + 3*F[n-1] + F[n])
       end
       integ = integ * h
   end
   return integ
end

"""
    integ1(n::Int64,h::Float64,F)

    integ1 performs 1st order left boundary integration method for datavector `F` with `n` points with a distance of `h` between each of them.
    If n > 1 and F = 1, F will be interpreted as a constant function over the n points.
"""
function integ1(n::Int64,h::Float64,F)
    if n != length(F) && length(F) == 1
       F = F*ones(n)
   elseif n != length(F)
       println("n does not math number of data points")
       return NaN
    end
 
    if n <= 1
       return 0
    elseif n==2
       return h*F[1]
    else
        return h*(sum(F)-F[end])
    end
end

"""
    integ2(n::Int64,h::Float64,F)

    integ2 performs 2nd order integration (trapezoid) method for datavector `F` with `n` points with a distance of `h` between each of them.
    If n > 1 and F = 1, F will be interpreted as a constant function over the n points.
"""
function integ2(n::Int64,h::Float64,F)
    if n != length(F) && length(F) == 1
       F = F*ones(n)
   elseif n != length(F)
       println("n does not math number of data points")
       return NaN
    end
 
    if n <= 1
       return 0
    else
        return h*(sum(F) - 0.5*F[1] - 0.5*F[end])
    end
end


"""
    CumInteg4(n::Int64,h::Float64,F)

    CumInteg4 calculates the cumulative integrals ``\\int_0^t F(t)dt \\forall t\\in(0,T)`` of `F` using a 4th order integration method.
    `n` = number of datapoints,
    `h` = distance between datapoints.
"""
function CumInteg4(n::Int64,h::Float64,F)
   if n != length(F) && length(F) == 1
      # F is assumed to be constant
      F = F*ones(n)
   elseif n > length(F)
       println("n does not math number of data points (CumInteg4)")
       println("n = ",n)
       println("F = ",F)
       sleep(5)
       return NaN
   end
   CumInt = zeros(n)
   if n < 1
   elseif n==1
       CumInt[1] = 0
   elseif n==2
       # trapezoidal rule
       CumInt[1] = 0
       CumInt[2] = h*0.5*(F[1] + F[2])
   else
       CumInt[1] = 0.0
       CumInt[2] = 0.5*(F[1] + F[2])
       CumInt[3] = 1/3 * (F[1] + 4*F[2] + F[3])

       for ii = 5:2:n
           CumInt[ii-1] = CumInt[ii-4] + 3/8 * (F[ii-4] + 3*F[ii-3] + 3*F[ii-2] + F[ii-1])
           CumInt[ii] = CumInt[ii-2] + 1/3*(F[ii-2] + 4*F[ii-1] + F[ii])
       end
       if iseven(n)
           CumInt[n] = CumInt[n-3] + 3/8 * (F[n-3] + 3*F[n-2] + 3*F[n-1] + F[n])
       end
       CumInt = CumInt * h
   end
   return CumInt
end


"""
    interp1(xpt,ypt;method = "Linear",extrapolation = Flat())

    One dimensional interpolation using either linear interpolation or BSplines
"""
function interp1(xpt,ypt;method = "Linear",extrapolation = Flat())
    if length(ypt) <=2
        method = "Linear"
    end
    if method == "Linear"
        f = interpolate((xpt,),ypt,Gridded(Linear()))
        f = extrapolate(f,extrapolation)
    elseif method == "Spline"
        f1 = interpolate(ypt, BSpline(Cubic(Natural(OnCell()))))
        f1 = extrapolate(f1,extrapolation)

        x0 = xpt[1]
        xdelta = xpt[2]-xpt[1]

        f = t-> f1(1 .+ (t.-x0)./xdelta )
    end
    return f
end

"""
    interp2(xpt,ypt,zpt)

    Two dimensional interpolation using BSplines
"""
function interp2(xpt,ypt,zpt)
    f1 = interpolate(zpt, BSpline(Cubic(Natural(OnCell()))))
    f1 = extrapolate(f1,Line())

    x0 = xpt[1]
    xdelta = xpt[2]-xpt[1]
    y0 = ypt[1]
    ydelta = ypt[2]-ypt[1]

    f = (t,s) -> f1(1 .+ (t.-x0)./xdelta,1 .+ (s.-y0)./ydelta )
    return f
end


"""
    SavitskyGolaySmoothing(x::Vector, windowSize::Integer, polyOrder::Integer; deriv::Integer=0)

    Smooth a vector x using the Savitsky-Golay method
"""
function SavitskyGolaySmoothing(x::Vector, windowSize::Integer, polyOrder::Integer; deriv::Integer=0)

	#Some error checking
	@assert isodd(windowSize) "Window size must be an odd integer."
	@assert polyOrder < windowSize "Polynomial order must me less than window size."

	halfWindow = floor(Int,windowSize/2)

	#Setup the S matrix of basis vectors.
	S = zeros(windowSize, polyOrder+1)
	for ct = 0:polyOrder
		for w = 1:windowSize
			S[w,ct+1] = (w-1-halfWindow)^(ct)
		end
	end

	#Compute the filter coefficients for all orders
	#From the scipy code it seems pinv(S) and taking rows should be enough
	G = S*pinv(S'*S)

	#Slice out the derivative order we want
	filterCoeffs = G[:,deriv+1] * factorial(deriv);

	#Pad the signal with the endpoints and convolve with filter
	paddedX = vcat(x[1]*ones(halfWindow), x, x[end]*ones(halfWindow))
	y = conv(filterCoeffs[end:-1:1], paddedX)

	#Return the valid midsection
	return y[2*halfWindow+1:end-2*halfWindow]
end


"""
    ParaAdjust(h_global::Float64,Para::Dict) --> VOID function

    ParaAdjust adjust all relevant elements of the Parameter dictionary `Para` to the step-size `h_global`.
    No Output as entries of `Para` get overwritten
"""
function ParaAdjust(h_global::Float64,Para::Dict)
   Para["hstep"] = h_global
   Para["nTime"] = floor(Int,Para["T"]/Para["hstep"]) + 1           # number of breakpoints on the time axis
   Para["T"] = Para["hstep"] * (Para["nTime"]-1)
   Para["nAge"] = floor(Int,Para["T"]/Para["hstep"]) + 1        # number of breakpoints on the age axis
   Para["tmesh"] = collect(0:Para["hstep"]:Para["T"])
   Para["amesh"] = collect(0:Para["hstep"]:Para["T"])

   Para["OptiIter"] = 1
end



"""
    LoadVariables(Para,Results)

    Load inital profiles for the control variables from the Results dictionary and use interpolation if necessary
"""
function LoadVariables(Para,Results)

    Con, Con_dist = ConInterpol(Results["Con"],Results["Con_dist"],Para)

    for kk in Para["Reduce1"]
        Con[:,:,kk] .= 0
    end
    for kk in Para["Reduce2"]
        Con_dist[:,:,kk] .= 0
    end

    return Con, Con_dist
end

"""
    ConInterpol(Con_dist::Array{Float64,3},Con::Array{Float64,3},Para::Dict)

    Given Controls `Con_dist` and `Con`, ConInterpol performs an interpolation of these fields to fit the passed Parameter dict `Para`.

    Outputs: New Controls 'Con_dist_new' and `Con_new`
"""
function ConInterpol(Con::Array{Float64,3},Con_dist::Array{Float64,3},Para::Dict)

    nAge_old, nTime_old, nCon_dist = size(Con_dist)
    nCon = size(Con,3)

    tmesh_old = collect(0:(1/(nTime_old-1)):1) * Para["T"]
    amesh_old = collect(0:(1/(nAge_old-1)):1)  * Para["T"]

    Con_dist_new = zeros(Para["nAge"],Para["nTime"],nCon_dist)
    Con_new = zeros(1,Para["nTime"],nCon)

    #----------------------------
    #   Interpolation
    #----------------------------

    if nTime_old < Para["nTime"]
        # Interpolation for distributed controls
        # Interpolation along the characteristic lines
        factor = floor(Int,(Para["nTime"]-1)/(nTime_old-1))       # (1/(nTime_old-1))/(1/(Para["nTime-1))
        for jj = 1:nCon_dist
            #   Along existing characteristics
            for kk = 1:(nAge_old-1)
                kk2 = 1+factor*(kk-1)
                f1 = interp1(tmesh_old[kk:end],Con_dist[kk,kk:end,jj],method = "Spline",extrapolation = Line())
                Con_dist_new[kk2,kk2:end,jj] = f1(Para["tmesh"][kk2:end])
            end
            Con_dist_new[end,end,jj] = Con_dist[end,end,jj]

            for kk = 1:Para["nTime"]-2
                f1 = interp1(Para["tmesh"][kk:2:end],diag(Con_dist_new[:,:,jj],kk-1)[1:2:end],method = "Spline",extrapolation = Line())
                for mm = 1:Para["nTime"]+1-kk
                    Con_dist_new[mm,mm+kk-1,jj] = f1(Para["tmesh"][mm+kk-1])
                end
            end
            
            f1 = interp1(tmesh_old,Con_dist[:,end,jj],method = "Spline",extrapolation = Line())
            Con_dist_new[:,end,jj] = f1(Para["tmesh"])
        end

        # Interpolation for non-distributed controls
        for jj = 1:nCon
            f1 = interp1(tmesh_old,Con[1,:,jj],method = "Spline",extrapolation = Line())
            Con_new[1,:,jj] = f1(Para["tmesh"])
        end

    elseif nTime_old > Para["nTime"]
        for kk = 1:nCon_dist
            Con_dist_new[:,:,kk] = evalgrid(Spline2D(amesh_old,tmesh_old,Con_dist[:,:,kk];kx=2,ky=2),Para["amesh"],Para["tmesh"])
        end
        for kk = 1:nCon
            Con_new[1,:,kk] =evaluate(Spline1D(tmesh_old,Con[1,:,kk],k=2),Para["tmesh"])
        end
    else
        return Con, Con_dist
    end


    #------------------------------
    #   Control Mapping
    #------------------------------
    ConMapping(Con_new,Para)
    for ii = 1:Para["nAge"]
        ConMapping_dist(Con_dist_new,ii,Para)
    end

    #----------------------------
    #   Smoothing
    #----------------------------
    if Para["ConSmooth"] == true
        ConSmooth(Con_new,Con_dist_new,Para)
        ConMapping(Con_new,Para)
        for ii = 1:Para["nAge"]
            ConMapping_dist(Con_dist_new,ii,Para)
        end
    end

    return Con_new, Con_dist_new
end

"""
    ConSmooth(Con::Array{Float64,3},Con_dist::Array{Float64,3},Para::Dict)

    Smooth the profiles of the control variables using a SavitskyGolaySmoother
"""
function ConSmooth(Con::Array{Float64,3},Con_dist::Array{Float64,3},Para::Dict)

    #---------------------------------------------------------------------------
    # Smoothing with SavitskyGolaySmoothing
    # Choose window with length of 10 years
    window = floor(Int,(Para["nTime"]-1)/Para["T"]*10)+1
    if window%2 == 0
        window = window-1
    end

    # Smoothing along vintages
    for kk = 1:Para["nCon"]
        Con[1,:,kk] = SavitskyGolaySmoothing(Con[1,:,kk],window,min(3,window-1))
    end

    for kk = 1:Para["nCon_dist"]
        for jj = 1:(Para["nTime"]-4)
            if window <= Para["nTime"] - jj +1
                window2 = window
            else
                window2 = Para["nTime"]-jj+1
                if window2%2 == 0
                    window2 = window2 - 1
                end
            end
            Con_dist[jj,jj:end,kk] = SavitskyGolaySmoothing(Con_dist[jj,jj:end,kk],window2,min(3,window2-1))
        end
    end
end

"""
    GradSmooth(dHam::Array{Float64,3},dHam_dist::Array{Float64,3},Para::Dict)

    Smooth the gradient by using the median value of a moving window.
"""
function GradSmooth(dHam::Array{Float64,3},dHam_dist::Array{Float64,3},Para::Dict)
    for ii = 3:(-1):1
        maxLag = ii
        for kk = 1:Para["nCon_dist"]
            for ii = 1:Para["nTime"]
                for jj = (ii+maxLag):(Para["nTime"]-maxLag)
                    first = max(ii,jj-maxLag)
                    last = min(Para["nTime"],jj+maxLag)
                    dHam_dist[ii,jj,kk] = quantile(dHam_dist[ii,first:last,kk],0.5)
                end
            end
        end
        for kk = 1:Para["nCon"]
            for jj = (1+maxLag):(Para["nTime"]-maxLag)
                first = max(1,jj-maxLag)
                last = min(Para["nTime"],jj+maxLag)
                dHam[1,jj,kk] = quantile(dHam[1,first:last,kk],0.5)
            end
        end
    end
end

function printSuccess(Para::Dict,Obj::Float64,Step)
    @printf("%g Success: ObjVal = %.4f, OptiStep=%.8f, OptiIter = %g, hStep = %.3f \n",Para["LineIter"], Obj, Step, Para["OptiIter"], Para["hstep"])
end

function printNoSuccess(Para::Dict,Obj::Float64,Step)
    @printf("%g No success: ObjVal = %.4f, OptiStep=%.8f, OptiIter = %g, hStep = %.3f \n",Para["LineIter"],Obj, Step, Para["OptiIter"], Para["hstep"])
end
