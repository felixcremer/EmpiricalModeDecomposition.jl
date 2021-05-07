using Printf

"""
Iterator for the sifting algorithm.
"""
mutable struct SiftIterable{solT <: AbstractVector, vecT <: AbstractVector, numT <: Real}
    "Measurements/IMF"
    y::solT
    "Time series values"
    x::vecT
    "Positions of the local maxima"
    maxs::Vector{Int}
    "Positions of the local minima"
    mins::Vector{Int}
    "Indices of the zero crossings"
    crosses::Vector{Int}
    "Tolerance"
    tol::numT
    "Residuals used to check stopping criteria"
    residual::numT
    "Maximum number of iteration"
    maxiter::Int
    "Number of steps of a stable sift after which the sifting is aborted"
    stop_steps::Int
    "Number of steps on which the number of zero crossings was fix"
    fix_steps::Int
end

@inline converged(it::SiftIterable) = it.residual ≤ var(it.y) * it.tol

@inline start(it::SiftIterable) = 0

@inline done(it::SiftIterable, iteration::Int) = iteration ≥ it.maxiter || converged(it)

# kernel sifting algorithm
function iterate(it::SiftIterable, iteration::Int=start(it))
    # Check for termination first
    if done(it, iteration)
        return nothing
    end

    localmaxmin!(it.y, it.maxs, it.mins)
    maxlen = length(it.maxs)
    minlen = length(it.mins)

    it.fix_steps == it.stop_steps && return nothing

    zerocrossing!(it.y, it.crosses)

    abs(length(it.crosses) - maxlen - minlen) <= 1 && (it.fix_steps += 1)
    if maxlen < 4 || minlen < 4
        return nothing
    end
    smin = get_edgepoint(it.y, it.x, it.mins, first, isless)
    smax = get_edgepoint(it.y, it.x, it.maxs, first, !isless)
    emin = get_edgepoint(it.y, it.x, it.mins, last, isless)
    emax = get_edgepoint(it.y, it.x, it.maxs, last, !isless)

    # TODO: Format for clarity?
    maxTS = interpolate([first(it.x); it.x[it.maxs]; last(it.x)], [smax; it.y[it.maxs]; emax], it.x, DierckXInterp())
    minTS = interpolate([first(it.x); it.x[it.mins]; last(it.x)], [smin; it.y[it.mins]; emin], it.x, DierckXInterp())
    subs = 0.5*(maxTS + minTS)
    it.residual = sum(abs, subs)

    it.y .-= subs

    # Return the residual at item and iteration number as state
    it.residual, iteration + 1
end

# Utility functions

function sift_iterator(y, x;
                      reltol::Real = 0.1,
                      maxiter::Int = 6,
                      stop_steps::Int = 4)
    maxs, mins = Int[], Int[]
    crosses = Int[]
    residual = sum(abs, y)
    tolerance = 0.1
    fix_steps = 0
    # Return the iterable
    return SiftIterable(y, x, maxs, mins, crosses,
        tolerance, residual, maxiter, stop_steps, fix_steps)
end

"""
    sift(y, x=1:length(y); kwargs...) -> imf

Same as [`sift!`](@ref), but allocates a solution vector `imf` initialized with zeros.
"""
sift(y, x=1:length(y); kwargs...) = sift!(zeros(eltype(y), size(y)), y, x; kwargs...)

"""
    sift!(imf, y, x; kwargs...) -> imf

# Arguments

- `imf`: intrinsic mode function, will be updated in-place;
- `y`: data measurements;
- `x`: time stamp series.

## Keywords

- `reltol::Real = 0.1`: relative tolerance for the stopping condition;
- `maxiter::Int = 10`: maximum number of iterations;
- `verbose::Bool = false`: print method information;
- `log::Bool = false`: keep track of the residual in each iteration.

# Output

**if `log` is `false`**

- `imf`: intrinsic mode function.
"""
function sift!(imf, y, x;
             reltol::Real = 0.1,
             maxiter::Int = 6,
             stop_steps::Int = 4,
             log::Bool = false,
             verbose::Bool = false,
             kwargs...)

    imf .= y
    # Actually perform sifting
    iterable = sift_iterator(imf, x; reltol, maxiter, stop_steps)

    for (iteration, item) = enumerate(iterable)
        verbose && @printf("%3d\t%1.2e\n", iteration, iterable.residual)
    end

    verbose && println()

    iterable.y
end