using TimeseriesSurrogates
using Plots
# Idea is, that I can set a function and a number of bootstraps with surrogate method
# And the bootstrap function (or macro) redos the function with the result replaced by the surrogate of the result.


function bootstrap(ts, surrfun, fun, ntimes=10)
    org = fun(ts)
    @show size(org)
    surr_start = ts .- org

    surrogates = typeof(org)[]
    for i in  1:ntimes
        surr = surrfun(surr_start)

        push!(surrogates, fun(surr + org))
    end
    return org, surrogates
end


function plot_surr(dates, org, surrs)
    pyplot()
    surr25 = zero(org)
    surr75 = zero(org)
    surrmed = zero(org)
    for i in eachindex(org)
       surr75[i] = quantile([surr[i] for surr in surrs], .75)
       surr25[i] = quantile([surr[i] for surr in surrs], .25)
       surrmed[i] = median([surr[i] for surr in surrs])
    end
    plot(dates, surr25, fill=(surr75,:orange), line=nothing)
    plot!(dates, org, color=:blue)
    plot!(dates, surrmed, color=:red)
    plot!(dates, surr25)
    plot!(dates, surr75)

    #return surr25, surr75
end

macro bootstrap(nstraps, surrogate_fun, ex)
    :(for i=1:$(esc(nstraps))
        $(esc(ex))
    end)
end
export bootstrap

deref(x) = x
deref(x::Array{T,0}) where T = x[]
macro bootstrap(ntimes, surr, expr)
    return quote
        org = $(esc(expr))
        @show expr
        @show org
        x = [($(esc(expr))) for i = 1:$(esc(ntimes))]
        @show x
        mu = mean(x)
        sigma = std(x)
        s = -log10.((/).(sigma, abs.(mu)))
        zip(mu, s) |> collect |> deref
    end
end
