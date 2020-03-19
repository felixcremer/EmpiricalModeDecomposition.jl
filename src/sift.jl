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
    @debug "Absolute sum of the data $s"
    crosses = Int[]

    state = SiftState(iter.yvec, iter.xvec,maxes, mins, crosses, s,0)
    return state ,state
end



function iterate(iter::SiftIterable, state::SiftState)
    localmaxmin!(state.yvec, state.maxes, state.mins)
    maxlen = length(state.maxes)
    minlen = length(state.mins)
    @debug "The number of minima: $minlen"
    @debug "The number of maxima: $maxlen"
    state.fix_steps == iter.stop_steps && return nothing
    zerocrossing!(state.yvec,state.crosses)
    abs(length(state.crosses) - maxlen - minlen) <=1 && (state.fix_steps +=1)
    if maxlen<1 || minlen<1
        return nothing
    end
    @debug state.yvec, state.xvec, state.mins, state.maxes
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
    @debug [first(state.xvec); state.xvec[state.maxes]; last(state.xvec)]
    @debug [smax; state.yvec[state.maxes]; emax]
    maxTS = EmpiricalModeDecomposition.interpolate([first(state.xvec); state.xvec[state.maxes]; last(state.xvec)],[smax; state.yvec[state.maxes]; emax],state.xvec,DataInterp())
    #state.yvec[1] = smin
    #state.yvec[end] = emin
    minTS = EmpiricalModeDecomposition.interpolate([first(state.xvec); state.xvec[state.mins]; last(state.xvec)],[smin; state.yvec[state.mins]; emin],state.xvec,DataInterp())
    subs = (maxTS+minTS)/2
    state.s =sum(abs,subs)
    @debug "Absolute sum of the not yet decomposed data $(state.s)"
    state.yvec = state.yvec - subs
    return state, state
end


"""
    sift(y, xvec)

    Sift the vector y whose points have x coordinates given by xvec.

"""
function sift(yvec, xvec=1:length(yvec), tol=0.1)
    ϵ = var(yvec) * tol
    @debug ϵ
    stop(state) = state.s <= ϵ
    imf=nothing
    num_steps = 0
    for (i, step) in enumerate(halt(SiftIterable(yvec, xvec, 4), stop))
        @debug sum(abs, step.yvec)
        imf= step.yvec
        num_steps = i
        @debug num_steps, step.s
    end
    @debug num_steps
    imf
end