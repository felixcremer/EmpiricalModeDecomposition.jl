"""
    localmaxmin!(x, maxes, mins)

Detect the local extrema of x.
Push the maxima into maxes and the minima into mins.
"""
function localmaxmin!(y, maxes::Vector{Int}, mins::Vector{Int})
    empty!(maxes)
    empty!(mins)
    for i in 2:length(y)-1
        if y[i+1] < y[i] > y[i-1]
            push!(maxes,i)
        elseif y[i+1] > y[i] < y[i-1]
            push!(mins, i)
        end
    end
end


"""
    zerocrossing!(y, crosses)

Compute the indices of zero crossings of vector `y`.
It searches for elements which are either zero or near a signflip and pushes the
indices into crosses.
"""
function zerocrossing!(y, crosses)
    empty!(crosses)
    for i ∈ 1:(length(y)-1)
        if y[i] == zero(y[i]) || sign(y[i]) * sign(y[i+1]) == -1
            push!(crosses, i)
        end
    end
end

"""
    get_edgepoint(y, xvec, extremas, pos, comp)

Compute the edgepoint which should be used as the extrema on the edge for the
spline computation.
"""
function get_edgepoint(y, xvec, extremas, pos, comp)
    if pos == first
        index = [1,2]
    elseif pos == last
        index = [length(extremas) - 1, length(extremas)]
    end
    # the x values must be embedded into a tuple
    knots = (xvec[extremas[index]],)
    itp = Interpolations.interpolate(knots, y[extremas[index]], Gridded(Linear()))
    expf = extrapolate(itp, Line())
    edgepoint = expf(pos(xvec))
    if comp(edgepoint, pos(y))
        edgepoint
    else
        edgepoint = pos(y)
    end
end

"""
    ismonotonic(x::AbstractVector)

Check whether `x` is monotonic.
"""
ismonotonic(x::AbstractVector) = issorted(x) || issorted(x, rev=true)

abstract type InterpMethod end
struct DierckXInterp <: InterpMethod end
#immutable InterpolationsInterp <: InterpMethod end

function interpolate(knotxvals::Vector, knotyvals::Vector, predictxvals::AbstractVector, m::DierckXInterp, k=3)
    spl = Dierckx.Spline1D(knotxvals, knotyvals, k=k)
    #@show spl, predictxvals
    Dierckx.evaluate(spl,predictxvals)
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
