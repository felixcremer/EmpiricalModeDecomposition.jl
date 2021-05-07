module EmpiricalModeDecomposition

using Interpolations
using Dierckx
using IterTools
using Random
using Statistics
using Base.Iterators

import Base: iterate, IteratorSize, eltype

export EMDIterable
export sift, sift!
export emd, eemd, ceemd, maketestdata

include("testdata.jl")
include("utils.jl")
include("sift.jl")


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
    if imf_prev[2]
        return nothing
    elseif sum(abs, imf_prev[1]) > 0.0 && !ismonotonic(imf_prev[1])
        imf = sift(imf_prev[1], iter.xvec; reltol=0.1)
        return imf, (imf_prev[1]-imf, false)
    else
        return imf_prev[1], (imf_prev[1], true)
    end
end

Base.IteratorSize(::Type{EMDIterable{U,V}}) where {U,V} = Base.SizeUnknown()


"""
    emd(measurements, xvec, numimfs=6)

Return the Intrinsic Mode Functions (IMF) and the residual of the Empirical Mode
Decomposition (EMD) of the `measurements` given on time steps given in `xvec`.
"""
function emd(measurements, xvec, num_imfs=6)
    collect(take(EMDIterable(measurements, xvec), num_imfs))
end



"""
    eemd(measurements, xvec, numtrails=100)

Return the Intrinsic Mode Functions and the residual of the ensemble Empirical Mode
Decomposition (EMD) of the `measurements` given on time steps `xvec`.
"""
function eemd(measurements, xvec, numtrails=100, num_imfs=6)
    random = randn(length(measurements))
    imfs_mean = emd(measurements .+ random, xvec, num_imfs)
    num_imfs = length(imfs_mean)
    for i in 1:numtrails
        randn!(random)
        imfs = emd(measurements .+ random, xvec, num_imfs)
        if length(imfs) < length(imfs_mean)
            imfs_mean[length(imfs)] = sum(imfs_mean[length(imfs):end])
            imfs_mean = imfs_mean[1:length(imfs)]
        elseif length(imfs) > length(imfs_mean)
            imfs[length(imfs_mean)] = sum(imfs[length(imfs_mean):end])
            imfs = imfs[1:length(imfs_mean)]
        end
        imfs_mean .+= imfs
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
    imf0 = mean([sift(iter.yvec+noise, iter.xvec; reltol=0.1) for noise in iter.noise_ens])
    iter_ens = EMDIterable.(iter.noise_ens, [iter.xvec])
    state_ens = Tuple[]
    imf_state_ens = iterate.(iter_ens)
    imf0, CEEMDState(iter.yvec-imf0, iter_ens, imf_state_ens, false)
end

function Base.iterate(iter::CEEMDIterable, state::CEEMDState)

    vstop = var(iter.yvec)*1e-10

    if state.finished

        return nothing

    elseif sum(abs,state.yvec)>vstop && !ismonotonic(state.yvec)

        # TODO: Format for clarity?
        imf = vec(median(hcat([sift(state.yvec+noise[1], iter.xvec; reltol=0.1) for noise in state.imf_state_ens]...),dims = 2))

        for iens in 1:length(state.iter_ens)
            r = iterate(state.iter_ens[iens], state.imf_state_ens[iens][2])
            if isnothing(r)
                fill!(state.imf_state_ens[iens][1], 0)
            else
                imf_noise, ensstate = r
                state.imf_state_ens[iens] = imf_noise, ensstate
            end
        end

        newstate = CEEMDState(state.yvec-imf, state.iter_ens, state.imf_state_ens, false)
        return imf,newstate
    else
        @show state.yvec
        return state.yvec,CEEMDState(state.yvec, state.iter_ens, state.imf_state_ens, true)
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
    maxs = Int[]
    mins = Int[]
    localmaxmin!(imf, maxs, mins)
    interpolate(xs[maxs], imf[maxs], xs, DierckXInterp())
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
