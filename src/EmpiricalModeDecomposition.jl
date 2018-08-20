module EmpiricalModeDecomposition
using Interpolations
using Dierckx
using IterTools
using Random


export emd, eemd, ceemd, maketestdata

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
    #@show startline
    startline<first(y) ? startline : first(y)
end

function startmax(y,xvec, maxes)
    startline = EmpiricalModeDecomposition.interpolate(xvec[maxes[1:2]],y[maxes[1:2]],[1],DierckXInterp(),1)[1]
    #@show startline
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

import Base.iterate, Base.IteratorSize

export EMDIterable, SiftIterable

struct SiftIterable{T<:AbstractVector,U<:AbstractVector}
    yvec ::T
    xvec ::U
end

mutable struct SiftState
    yvec
    xvec
    maxes::Vector{Int}
    mins::Vector{Int}
    s
end

Base.IteratorSize(::Type{SiftIterable}) = Base.SizeUnknown()


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
    #@show ϵ
    stop(state) = state.s <= ϵ
    imf=nothing
    for step in halt(SiftIterable(yvec, xvec), stop)
        #@show sum(abs, step.yvec)
        imf= step.yvec
    end
    imf
end

function ceemd(measurements, xvec; num_imfs=6, numtrails=100, β=0.02, noise_ens = [β .* randn(length(xvec)) for i in 1:numtrails])
  CEEMDIterable(measurements,xvec,noise_ens)
end

struct EMDIterable{U<:AbstractVector,V<:AbstractVector}
    yvec::U
    xvec::V
end

Base.IteratorSize(::Type{EMDIterable{U,V}}) where {U,V} = Base.SizeUnknown()

import EmpiricalModeDecomposition: sift, ismonotonic
using Statistics
struct CEEMDIterable{U<:AbstractVector,V<:AbstractVector,T<:AbstractVector}
    yvec::U
    xvec::V
    noise_ens::T
end
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

    imf = mean([sift(state.yvec+noise[1],iter.xvec,0.1) for noise in state.imf_state_ens])

    for iens in 1:length(state.iter_ens)
      r = iterate(state.iter_ens[iens],state.imf_state_ens[iens][2])
      if r == nothing
        fill!(state.imf_state_ens[iens][1], 0)
      else
        imf, ensstate = r
        state.imf_state_ens[iens] = imf, ensstate
      end
    end

    newstate = CEEMDState(state.yvec-imf,state.iter_ens,state.imf_state_ens,false)
    return imf,newstate
  else
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


function maketestdata(seed)
  Random.seed!(seed)
  ## simulate data of length....
  N_tim = 240
  NpY   = 24         # samples/year
  t     = 1:N_tim
  t     = t./NpY # your time vector

  # constant seasonal cycle
  A   = 2          # amplitude
  phi = 13*pi/12   # initial phase
  S   = A*cos.(2*pi*t+phi)

  # generate a linear trend
  T = 0.1 + 0.0002.*t


  # some other oscillation
  a = 0.2
  b = 0.1
  C = (1+a*cos.(b*2*pi*t))

  # simple (coloured) noise
  φ = 0.3 # strengh of autocorrelation in noise
  E = randn(N_tim).*0.1

  for i = 2:N_tim
    E[i] = φ*E[i-1]+(1-φ)*E[i]
  end

  X = S.*C + 2.*E
  t,X,E
end

end # module
