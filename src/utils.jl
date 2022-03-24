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
    if pos == first
        if length(extremas) > 1
            index = [extremas[1], extremas[2]]
        else
            index = [1, extremas...]
        end
    elseif pos == last
        if length(extremas) > 1
            index = [extremas[end-1], extremas[end]]
        else
            index = [extremas..., length(xvec)]
        end
    end
    @debug index
    knots = (xvec[index],)
    itp = Interpolations.interpolate(knots,y[index], Gridded(Linear()))
    expf = Interpolations.extrapolate(itp, Interpolations.Line())
    edgepoint = expf(pos(xvec))
    @debug edgepoint
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
struct DataInterp <: InterpMethod end

function interpolate(knotxvals::Vector,knotyvals::Vector,predictxvals::AbstractVector,m::DierckXInterp, k=3)
    spl = Dierckx.Spline1D(knotxvals, knotyvals, k=k)
    @debug spl, predictxvals
    Dierckx.evaluate(spl,predictxvals)
end

function interpolate(knotxvals, knotyvals, predictxvals, m::DataInterp)
    spl=DataInterpolations.CubicSpline(knotyvals, knotxvals)
    @debug "Interpolation" spl
    return spl.(predictxvals)
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
