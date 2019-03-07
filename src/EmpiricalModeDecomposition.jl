module EmpiricalModeDecomposition
using Interpolations
using Dierckx
using IterTools
using Random
using Statistics
import Base.iterate, Base.IteratorSize

export EMDIterable, SiftIterable
export emd, eemd, ceemd, maketestdata

include("testdata.jl")

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


"""
zerocrossing!(y, crosses)

Compute the indices of zerocrossings of a vector
It searches for elements which are either zero or near a signflip
and pushes the indices into crosses.
"""
function zerocrossing!(y, crosses)
    empty!(crosses)
    for i ∈ 1:(length(y)-1)
        if y[i] == zero(y[i]) || sign(y[i]) * sign(y[i+1]) == -1
            push!(crosses, i)#
        end
    end
end

"""
get_edgepoint(y, xvec, extremas, pos, comp)

Compute the edgepoint which should be used as the extrema on the edge for the spline computation.
"""
function get_edgepoint(y, xvec, extremas, pos, comp)
 #the x values must be embedded into a tuple
    #@show knots
    if pos == first
        index = [1,2]
    elseif pos == last
        index = [length(extremas) - 1, length(extremas)]
    end
    knots = (xvec[extremas[index]],)
    itp = Interpolations.interpolate(knots,y[extremas[index]], Gridded(Linear()))
    expf = extrapolate(itp, Line())
    edgepoint = expf(pos(xvec))
    #@show edgepoint
    if comp(edgepoint, pos(y))
        edgepoint
    else
        edgepoint = pos(y)
    end
    return edgepoint
end

"""
ismonotonic(x::AbstractVector)

Check wether x is monotonic.
This means, every value is either larger or smaller than the preceding value.
"""
ismonotonic(x::AbstractVector{T}) where T = isfinite(foldl((x,y)->y>=x ? y : typemax(T), x, init=typemin(T))) || isfinite(foldl((x,y)->y<=x ? y : typemin(T), x, init=typemax(T)))

abstract type InterpMethod end
struct DierckXInterp <: InterpMethod end
#immutable InterpolationsInterp <: InterpMethod end

function interpolate(knotxvals::Vector,knotyvals::Vector,predictxvals::AbstractVector,m::DierckXInterp, k=3)
    spl = Dierckx.Spline1D(knotxvals, knotyvals, k=k)
    #@show spl, predictxvals
    Dierckx.evaluate(spl,predictxvals)
end


"""
SiftIterable{T, U}

Iterator for the sifting algorithm.
The time series values are an AbstractVector of type T,
and the time positions are an AbstractVector of type U.
Fields:
"""
struct SiftIterable{T<:AbstractVector,U<:AbstractVector}
    "Values of the time series"
    yvec ::T
    "Positions of the values"
    xvec ::U
    "Number of steps of a stable sift after which the sifting is aborted"
    stop_steps::Integer
end

"""
SiftState

Handle the  intermediate results of the sifting.
"""
mutable struct SiftState
    "Values of the time series"
    yvec
    "Position of the time steps"
    xvec
    "Positions of the local maxima"
    maxes::Vector{Int}
    "Positions of the local minima"
    mins::Vector{Int}
    "Indices of the zero crossings"
    crosses::Vector{Int}
    "sum of the abs value of yvec, this is used in the stopping criteria"
    s
    "Number of steps on which the number of zerocrossings was fix"
    fix_steps::Integer
end

Base.IteratorSize(::Type{SiftIterable}) = Base.SizeUnknown()


function iterate(iter::SiftIterable)
    maxes = Int[]
    mins = Int[]
    s =sum(abs, iter.yvec)
    crosses = Int[]

    state = SiftState(iter.yvec, iter.xvec,maxes, mins, crosses, s,0)
    return state ,state
end



function iterate(iter::SiftIterable, state::SiftState)
    localmaxmin!(state.yvec, state.maxes, state.mins)
    maxlen = length(state.maxes)
    minlen = length(state.mins)
    state.fix_steps == iter.stop_steps && return nothing
    zerocrossing!(state.yvec,state.crosses)
    abs(length(state.crosses) - maxlen - minlen) <=1 && (state.fix_steps +=1)
    if maxlen<4 || minlen<4
        return nothing
    end
    smin = get_edgepoint(state.yvec, state.xvec, state.mins, first, isless)
    smax = get_edgepoint(state.yvec, state.xvec, state.maxes, first, !isless)
    emin = get_edgepoint(state.yvec, state.xvec, state.mins, last, isless)
    emax = get_edgepoint(state.yvec, state.xvec, state.maxes, last, !isless)
    #state.yvec[1] = smax
    #state.yvec[end] =emax
    #pushfirst!(state.mins,1)
    #push!(state.mins, length(state.yvec))
    #pushfirst!(state.maxes,1)
    #push!(state.maxes, length(state.yvec))
    #@show [first(state.xvec); state.xvec[state.maxes]; last(state.xvec)]
    #@show [smax; state.yvec[state.maxes]; emax]
    maxTS = EmpiricalModeDecomposition.interpolate([first(state.xvec); state.xvec[state.maxes]; last(state.xvec)],[smax; state.yvec[state.maxes]; emax],state.xvec,DierckXInterp())
    #state.yvec[1] = smin
    #state.yvec[end] = emin
    minTS = EmpiricalModeDecomposition.interpolate([first(state.xvec); state.xvec[state.mins]; last(state.xvec)],[smin; state.yvec[state.mins]; emin],state.xvec,DierckXInterp())
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
    ϵ = var(yvec) * tol
    #@show ϵ
    stop(state) = state.s <= ϵ
    imf=nothing
    num_steps = 0
    for (i, step) in enumerate(halt(SiftIterable(yvec, xvec, 4), stop))
        #@show sum(abs, step.yvec)
        imf= step.yvec
        num_steps = i
    end
    #@show num_steps
    imf
end

"""
ceemd(measurements, xvec; num_imfs=6)

Compute the Complete Empirical Mode Decomposition
of the time series with values measurements and time steps xvec.
num_imfs is the Number of Intrinsic Mode Functions.
Returns a list of num_imfs + 1 Vectors of the same size as measurements.
"""
function ceemd(measurements, xvec; num_imfs=6, numtrails=100, β=0.04, noise_ens = [β*std(measurements) .* randn(length(xvec)) for i in 1:numtrails])
    imfs = collect(take(CEEMDIterable(measurements,xvec,noise_ens),num_imfs))
    @show size.(imfs)
    residual = measurements - sum(imfs)
    push!(imfs, residual)
    return imfs
end

"""
Iterator for the Empirical Mode Decomposition.
The time series values are an AbstractVector of type T,
and the time positions are an AbstractVector of type U.
"""
struct EMDIterable{T<:AbstractVector,U<:AbstractVector}
    yvec::T
    xvec::U
end

Base.IteratorSize(::Type{EMDIterable{U,V}}) where {U,V} = Base.SizeUnknown()

"""
Iterator for the Complete Empirical Mode Decomposition.
The time series values are an AbstractVector of type T,
and the time positions are an AbstractVector of type U.
"""
struct CEEMDIterable{U<:AbstractVector,V<:AbstractVector,T<:AbstractVector}
    "Vector with the values of the time series"
    yvec::U
    "Vector with the positions of the time series"
    xvec::V
    "Ensemble of noise which is a Vector of Vectors of the same size as yvec."
    noise_ens::T
end

"""
Intermediate results of the CEEMD Iteration
"""
struct CEEMDState
  yvec
  iter_ens
  imf_state_ens
  finished::Bool
end
Base.IteratorSize(::Type{CEEMDIterable{U,V,T}}) where {U,V,T} = Base.SizeUnknown()

function Base.iterate(iter::CEEMDIterable)
  imf0 = mean([sift(iter.yvec+noise,iter.xvec,0.1) for noise in iter.noise_ens])
  iter_ens = EMDIterable.(iter.noise_ens,[iter.xvec])
  state_ens = Tuple[]
  imf_state_ens = iterate.(iter_ens)
  imf0,CEEMDState(iter.yvec-imf0,iter_ens,imf_state_ens,false)
end

function Base.iterate(iter::CEEMDIterable,state::CEEMDState)

  vstop = var(iter.yvec)*1e-10

  if state.finished

    return nothing

  elseif sum(abs,state.yvec)>vstop && !ismonotonic(state.yvec)

      imf = vec(median(hcat([sift(state.yvec+noise[1],iter.xvec,0.1) for noise in state.imf_state_ens]...),dims = 2))

    for iens in 1:length(state.iter_ens)
      r = iterate(state.iter_ens[iens],state.imf_state_ens[iens][2])
      if r == nothing
        fill!(state.imf_state_ens[iens][1], 0)
      else
        imf_noise, ensstate = r
        state.imf_state_ens[iens] = imf_noise, ensstate
      end
    end

    newstate = CEEMDState(state.yvec-imf,state.iter_ens,state.imf_state_ens,false)
    return imf,newstate
  else
    @show state.yvec
    return state.yvec,CEEMDState(state.yvec,state.iter_ens,state.imf_state_ens,true)
  end
end



function iterate(iter::EMDIterable, imf_prev=(iter.yvec,false))
    if imf_prev[2]
      return nothing
    elseif sum(abs,imf_prev[1])>0.0 && !ismonotonic(imf_prev[1])
      imf = sift(imf_prev[1],iter.xvec, 0.1)
      return imf,(imf_prev[1]-imf,false)
    else
      return imf_prev[1],(imf_prev[1],true)
    end
end


using Base.Iterators
"""
emd(measurements, xvec)

Return the Intrinsic Mode Functions and
the residual of the Empirical Mode Decomposition of the measurements given on time steps given in xvec.
"""
function emd(measurements,xvec, num_imfs=6)
    collect(take(EMDIterable(measurements,xvec),num_imfs))
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
        #@show length(imfs_mean)

        randn!(random)
        imfs = EmpiricalModeDecomposition.emd(measurements .+ random, xvec, num_imfs)
        if length(imfs) < length(imfs_mean)
            imfs_mean[length(imfs)] = sum(imfs_mean[length(imfs):end])
            imfs_mean = imfs_mean[1:length(imfs)]
        elseif length(imfs) > length(imfs_mean)
            imfs[length(imfs_mean)] = sum(imfs[length(imfs_mean):end])
            imfs = imfs[1:length(imfs_mean)]
        end
        #@show length(imfs), length(imfs_mean)
        imfs_mean .+= imfs
    end

    imfs_mean ./= numtrails
end



function tone_masking(ys, xs, tone)
    ys_plus = ys .+ tone
    ys_minus = ys .- tone
    phi_plus = sift(ys_plus)
    phi_minus = sift(ys_minus)
    return (phi_plus .+ phi_minus) ./ 2
end


function iaestimation(imf, xs)
    r = abs.(imf)
    maxes = Int[]
    mins = Int[]
    localmaxmin!(imf, maxes, mins)
    EmpiricalModeDecomposition.interpolate(xs[maxes], imf[maxes], xs, DierckXInterp())
end


function iterAMremoval(imf, xvec)
    g = imf
    n=1
    b = zero(imf)
    while any(b .!=1) && n<=3
        b = iaestimation(imf, xvec)
        g = g ./ b
        n+=1
    end
    g
end




#
# function ceemd(measurements, xvec, num_imfs=6, numtrails=100, β=0.02, noise_ens = [β .* randn(length(xvec)) for i in 1:numtrails])
#     imfs = typeof(measurements)[]
#     #noise_imfs = [EMD(noise) for noise in noise_ens]
#     ycur = sift(measurements)
#     for i in 1:numtrails
#         ycur += sift(measurements .+ noise_ens[i],xvec)
#     end
#     ycur ./= numtrails
#     k=0
#     while length(imfs)<num_imfs && sum(abs, ycur) > 0.0 && !ismonotonic(ycur)
#         println(sum(ycur))
#         global k+=1
# #        plot(ycur)
#         y = sift(ycur, )
#         y = @parallel (+) for i in 1:numtrails
#             summand = ycur + 0.02*EmpiricalModeDecomposition.emd(randn(length(ycur)), xvec, k)[end]
#             EmpiricalModeDecomposition.sift(summand, xvec)
#         end
#         y ./=numtrails
#         println(maximum(y))
#         ycur = ycur-y
#         push!(imfs, ycur)
#     end
#     imfs
# end

end # module
