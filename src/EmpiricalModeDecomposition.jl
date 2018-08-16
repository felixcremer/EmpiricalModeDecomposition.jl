module EmpiricalModeDecomposition
using Interpolations
using Dierckx
using IterTools
using Random


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
    startline = interpolate(xvec[mins[1:2]],y[mins[1:2]],[1],DierckXInterp(),1)[1]
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

function endmin(y,xvec, maxes)
    startline = EmpiricalModeDecomposition.interpolate(xvec[maxes[end-1:end]],y[maxes[end-1:end]],[1],DierckXInterp(),1)[1]
    #@show startline
    startline<y[end] ? startline : y[end]
end

ismonotonic(x::AbstractArray{T}) where T = isfinite(foldl((x,y)->y>=x ? y : typemax(T), x, init=typemin(T))) || isfinite(foldl((x,y)->y<=x ? y : typemin(T), x, init=typemax(T)))

abstract type InterpMethod end
struct DierckXInterp <: InterpMethod end
#immutable InterpolationsInterp <: InterpMethod end

function interpolate(knotxvals::Vector,knotyvals::Vector,predictxvals::AbstractVector,m::DierckXInterp, k=3)
    spl = Dierckx.Spline1D(knotxvals, knotyvals, k=k)
    Dierckx.evaluate(spl,predictxvals)
end

import Base.iterate


struct SiftIterable
    yvec
    xvec
end

mutable struct SiftState
    yvec
    xvec
    maxes::Vector{Int}
    mins::Vector{Int}
    s
end



function iterate(iter::SiftIterable)
    maxes = Int[]
    mins = Int[]
    s =sum(abs, iter.yvec)
    state = SiftState(iter.yvec, iter.xvec,maxes, mins, s)
    return state,state
end



function iterate(iter::SiftIterable, state::SiftState)
    localmaxmin!(state.yvec, state.maxes, state.mins)
    if length(state.maxes)<4 || length(state.mins)<4
        return nothing
    end
    maxTS = EmpiricalModeDecomposition.interpolate(state.xvec[state.maxes],state.yvec[state.maxes],state.xvec,DierckXInterp())
    minTS = EmpiricalModeDecomposition.interpolate(state.xvec[state.mins],state.yvec[state.mins],state.xvec,DierckXInterp())
    subs = (maxTS+minTS)/2
    state.s =sum(abs,subs)
    state.yvec = state.yvec - subs
    return state, state
end


struct HaltingIterable{I, F}
    iter::I
    fun::F
end

function iterate(iter::HaltingIterable)
    next = iterate(iter.iter)
    return dispatch(iter, next)
end

function iterate(iter::HaltingIterable, (instruction, state))
    if instruction == :halt return nothing end
    next = iterate(iter.iter, state)
    return dispatch(iter, next)
end

function dispatch(iter::HaltingIterable, next)
    if next === nothing return nothing end
    return next[1], (iter.fun(next[1]) ? :halt : :continue, next[2])
end

halt(iter::I, fun::F) where {I,F} = HaltingIterable{I,F}(iter, fun)

"""
    sift(y, xvec)

    Sift the vector y whose points have x coordinates given by xvec.
"""
function sift(yvec, xvec=1:length(yvec), tol=0.1)
    ϵ = sum(abs, yvec) * tol
    @show ϵ
    stop(state) = state.s <= ϵ
    imf=nothing
    for step in halt(SiftIterable(yvec, xvec), stop)
        @show sum(abs, step.yvec)
        imf= step.yvec
    end
    imf
end



struct EMDIterable
    yvec
    xvec
end



function iterate(iter::EMDIterable, imf_prev=iter.yvec)
    imf = sift(imf_prev,state.xvec, 0.1)
    return imf
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
    while length(imfs)<num_imfs && sum(abs,ycur)>0.0 && !ismonotonic(ycur)
    #    println("While: ", length(imfs))
        #@show ycur
        y = sift(ycur,xvec)
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
function eemd(measurements, xvec, numtrails=100, num_imfs=6)
    random = randn(length(measurements))
    imfs_mean = emd(measurements .+ random, xvec, num_imfs)
    num_imfs = length(imfs_mean)
    for i in 1:numtrails
        @show length(imfs_mean)

        randn!(random)
        imfs = EmpiricalModeDecomposition.emd(measurements .+ random, xvec, num_imfs)
        if length(imfs) < length(imfs_mean)
            imfs_mean[length(imfs)] = sum(imfs_mean[length(imfs):end])
            imfs_mean = imfs_mean[1:length(imfs)]
        elseif length(imfs) > length(imfs_mean)
            imfs[length(imfs_mean)] = sum(imfs[length(imfs_mean):end])
            imfs = imfs[1:length(imfs_mean)]
        end
        @show length(imfs), length(imfs_mean)
        imfs_mean .+= imfs
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
