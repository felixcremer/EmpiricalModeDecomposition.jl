module EmpiricalModeDecomposition
using Dierckx

export emd, eemd, ceemd

"""
    localmaxmin!(x, maxes, mins)
Detect the local extrema of x.
Push the maxima into maxes and the minima into mins.
"""
function localmaxmin!(y,maxes::Vector{Int},mins::Vector{Int})
    empty!(maxes)
    empty!(mins)
    for i in 2:(length(y)-1)
        if y[i+1]<y[i]>y[i-1]
            push!(maxes,i)
        elseif y[i+1]>y[i]<y[i-1]
            push!(mins,i)
        end
    end
end

function startmin(y,xvec, mins)
    startline = EmpiricalModeDecomposition.interpolate(xvec[mins[1:2]],y[mins[1:2]],[1],DierckXInterp(),1)[1]
    @show startline
    startline<first(y) ? startline : first(y)
end

function startmax(y,xvec, maxes)
    startline = EmpiricalModeDecomposition.interpolate(xvec[maxes[1:2]],y[maxes[1:2]],[1],DierckXInterp(),1)[1]
    @show startline
    startline>first(y) ? startline : first(y)
end

function endmax(y,xvec, maxes)
    startline = EmpiricalModeDecomposition.interpolate(xvec[maxes[end-1:end]],y[maxes[end-1:end]],[1],DierckXInterp(),1)[1]
    #@show startline
    startline>y[end] ? startline : y[end]
end

function endmax(y,xvec, maxes)
    startline = EmpiricalModeDecomposition.interpolate(xvec[maxes[end-1:end]],y[maxes[end-1:end]],[1],DierckXInterp(),1)[1]
    #@show startline
    startline<y[end] ? startline : y[end]
end

ismonotonic{T}(x::AbstractArray{T})=isfinite(foldl((x,y)->y>=x ? y : typemax(T),typemin(T),x)) || isfinite(foldl((x,y)->y<=x ? y : typemin(T),typemax(T),x))

abstract type InterpMethod end
immutable DierckXInterp <: InterpMethod end
#immutable InterpolationsInterp <: InterpMethod end

function interpolate(knotxvals::Vector,knotyvals::Vector,predictxvals::AbstractVector,m::DierckXInterp, k=3)
    spl = Dierckx.Spline1D(knotxvals, knotyvals, k=k)
    Dierckx.evaluate(spl,predictxvals)
end

"""
    sift(y, xvec)

    Sift the vector y whose points have x coordinates given by xvec.
"""
function sift(y,xvec)
    #println("Sift")
    s = sum(abs,y)
    ϵ = s
    maxes = Int[]
    mins = Int[]
    while ϵ > 0.1*s
        #println("Sift While")
        EmpiricalModeDecomposition.localmaxmin!(y,maxes,mins)
        if length(maxes)<4 || length(mins)<4
            break
        end
        maxTS = EmpiricalModeDecomposition.interpolate(xvec[maxes],y[maxes],xvec,DierckXInterp())
        minTS = EmpiricalModeDecomposition.interpolate(xvec[mins],y[mins],xvec,DierckXInterp())
        subs = (maxTS+minTS)/2
        ϵ=sum(abs,subs)
        y = y - subs
    end
    y
end

"""
    emd(measurements, xvec)
Return the Intrinsic Mode Functions and
the residual of the Empirical Mode Decomposition of the measurements given on time steps given in xvec.
"""
function emd(measurements,xvec, num_imfs=100)
    imfs = typeof(measurements)[]
    ycur = measurements
    #println(length(ycur))
    #println(sum(abs,ycur)>0.0&& !EmpiricalModeDecomposition.ismonotonic(ycur))
    while length(imfs)<num_imfs && sum(abs,ycur)>0.0 && !EmpiricalModeDecomposition.ismonotonic(ycur)
    #    println("While: ", length(imfs))
        y = EmpiricalModeDecomposition.sift(ycur,xvec)
        push!(imfs,y)
        ycur = ycur-y
    end
    push!(imfs, ycur)
    imfs
end

"""
    eemd(measurements, xvec, numtrails=100)

    Return the Intrinsic Mode Functions and
    the residual of the ensemble Empirical Mode Decomposition of the measurements given on time steps xvec.
"""
function eemd(measurements, xvec, numtrails=100, num_imfs=5)
    imfs_mean  = @parallel (+) for i in 1:numtrails
                random = randn(length(xvec))
                imfs = EmpiricalModeDecomposition.emd(measurements+random, xvec, num_imfs)
                #println(length(imfs))
                imfs
            end

    imfs_mean ./= numtrails
end

function ceemd(measurements, xvec, num_imfs=6, numtrails=100, β=0.02)
    imfs = typeof(measurements)[]
    ycur = @parallel (+) for i in 1:numtrails
        EmpiricalModeDecomposition.sift(measurements+β*randn(length(xvec)),xvec)
    end
    k=0
    ycur ./= numtrails
    while length(imfs)<num_imfs && sum(abs, ycur) > 0.0 && !EmpiricalModeDecomposition.ismonotonic(ycur)
        println(sum(ycur))
        k+=1
#        plot(ycur)
        y = @parallel (+) for i in 1:numtrails
            summand = ycur + 0.02*EmpiricalModeDecomposition.emd(randn(length(ycur)), xvec, k)[end]
            EmpiricalModeDecomposition.sift(summand, xvec)
        end
        y ./=numtrails
        println(maximum(y))
        ycur = ycur-y
        push!(imfs, ycur)
    end
    imfs
end
end # module
