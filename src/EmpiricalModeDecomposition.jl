module EmpiricalModeDecomposition

using Interpolations
using Dierckx
using IterTools
using Random
using Statistics
using Base.Iterators
using DataInterpolations: DataInterpolations

import Base: iterate, IteratorSize, eltype

export EMDIterable, SiftIterable
export emd, eemd, ceemd, maketestdata

include("testdata.jl")
include("utils.jl")
include("sift.jl")
include("improvedceemdan.jl")


"""
Iterator for the Empirical Mode Decomposition.
The time series values are an AbstractVector of type `T`,
and the time positions are an AbstractVector of type `U`.
"""
struct EMDIterable{T<:AbstractVector, U<:AbstractVector}
    yvec::T
    xvec::U
end

eltype(::Type{EMDIterable{T,U}}) where {T,U} = T

function iterate(iter::EMDIterable, imf_prev=(iter.yvec,false))
    @debug sum(abs, imf_prev[1]), ismonotonic(imf_prev[1])

    if imf_prev[2]
      return nothing
    elseif sum(abs,imf_prev[1])>eps(eltype(iter.yvec))*length(iter.xvec) && !ismonotonic(imf_prev[1])
      imf = sift(imf_prev[1],iter.xvec, 0.1)
      @debug sum(abs, imf)
      @debug imf == imf_prev[1]
      return imf,(imf_prev[1]-imf,false)
    else
        @debug imf_prev
      return imf_prev[1],(imf_prev[1],true)
    end
end

Base.IteratorSize(::Type{EMDIterable{U,V}}) where {U,V} = Base.SizeUnknown()
eltype(::Type{EMDIterable{U,V}}) where {U,V} = U

"""
    emd(measurements, xvec, numimfs=6)

Return the Intrinsic Mode Functions (IMF) and the residual of the Empirical Mode
Decomposition (EMD) of the `measurements` given on time steps given in `xvec`.
"""
function emd(measurements,xvec, num_imfs=typemax(Int))
    imfs = collect(take(EMDIterable(measurements,xvec),num_imfs))
    #@show typeof(imfs)
    if !(sum(imfs) ≈ measurements)
      imfs[end] += measurements-sum(imfs)
    end
    imfs
  end



"""
    eemd(measurements, xvec, numtrails=100)

Return the Intrinsic Mode Functions and the residual of the ensemble Empirical Mode
Decomposition (EMD) of the `measurements` given on time steps `xvec`.
"""
function eemd(measurements, xvec; numtrails=100, num_imfs=ceil(Int, log(length(measurements))))
    random = randn(length(measurements))
    #@show num_imfs
    imfs_mean = [zero(measurements) for i in 1:num_imfs]
    #num_imfs = min(num_imfs, length(imfs_mean))
    #@show length(imfs_mean)
    for i in  1:numtrails
        @debug length(imfs_mean)
        randn!(random)
        imfs = EmpiricalModeDecomposition.emd(measurements .+ random, xvec, num_imfs)
        #if length(imfs) < length(imfs_mean)
        #    imfs_mean[length(imfs)] = sum(imfs_mean[length(imfs):end])
        #    imfs_mean = imfs_mean[1:length(imfs)]
        #elseif length(imfs) > length(imfs_mean)
        #    imfs[length(imfs_mean)] = sum(imfs[length(imfs_mean):end])
        #    imfs = imfs[1:length(imfs_mean)]
        #end
        #@show length(imfs), length(imfs_mean)
        imfs_mean[begin:length(imfs)] .+= imfs
    end

    imfs_mean ./= numtrails
end


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
  imf0 = mean([sift(iter.yvec+noise,iter.xvec) for noise in iter.noise_ens])
  #Why is iter.xvec wrapped in a vector here?
  iter_ens = EMDIterable.(iter.noise_ens,[iter.xvec])
  state_ens = Tuple[]
  imf_state_ens = iterate.(iter_ens)
  imf0,CEEMDState(iter.yvec-imf0,iter_ens,imf_state_ens,false)
end

function Base.iterate(iter::CEEMDIterable, state::CEEMDState)

    vstop = var(iter.yvec)*1e-10

    if state.finished

        return nothing

    elseif sum(abs,state.yvec)>vstop && !ismonotonic(state.yvec)

        # TODO: Format for clarity?
        imf = vec(median(hcat([sift(state.yvec+noise[1], iter.xvec, 0.1) for noise in state.imf_state_ens]...),dims = 2))

        for iens in 1:length(state.iter_ens)
            r = iterate(state.iter_ens[iens], state.imf_state_ens[iens][2])
            if isnothing(r)
                fill!(state.imf_state_ens[iens][1], 0)
            else
                imf_noise, ensstate = r
                state.imf_state_ens[iens] = imf_noise, ensstate
            end
        end

    newstate = CEEMDState(state.yvec-imf,state.iter_ens,state.imf_state_ens,false)
    return imf,newstate
  else
    @debug state.yvec
    return state.yvec,CEEMDState(state.yvec,state.iter_ens,state.imf_state_ens,true)
  end
end

"""
    ceemd(measurements, xvec; num_imfs=6)

Compute the Complete Empirical Mode Decomposition (CEMD) of the time series with values
`measurements` and time steps `xvec`.
`num_imfs` is the number of Intrinsic Mode Functions.
Returns a list of num_imfs + 1 vectors of the same size as measurements.
"""
function ceemd(measurements, xvec; num_imfs=6, numtrails=100, β=0.04,
    noise_ens = [β*std(measurements) .* randn(length(xvec)) for i in 1:numtrails])
    imfs = collect(take(CEEMDIterable(measurements, xvec, noise_ens), num_imfs))
    @show size.(imfs)
    residual = measurements - sum(imfs)
    push!(imfs, residual)
    return imfs
end

function tone_masking(ys, xs, tone)
    ys_plus = ys .+ tone
    ys_minus = ys .- tone
    phi_plus = sift(ys_plus)
    phi_minus = sift(ys_minus)
    @. 0.5 * (phi_plus + phi_minus)
end


function iaestimation(imf, xs)
    r = abs.(imf)
    maxes = Int[]
    mins = Int[]
    localmaxmin!(imf, maxes, mins)
    EmpiricalModeDecomposition.interpolate(xs[maxes], imf[maxes], xs, DataInterp())
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
