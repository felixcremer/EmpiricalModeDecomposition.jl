using TimeseriesSurrogates
using Plots
using LibEEMD
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

function surrmetrics(ts, x)
    surrmetrics(ts, x, y->ceemdan(y, 2)[:,end-1])
end

function surrmetrics(ts, x, fun)
    org, surr = bootstrap(ts, iaaft, fun,200)
    surr, surrmetrics(surr)
end

function surrmetrics(surrs)
    @show size(surrs[1])
    @show typeof(surrs[1])
    surr25 = zero(surrs[1])
    surr75 = zero(surrs[1])
    surr05 = zero(surrs[1])
    surr95 = zero(surrs[1])
    surrmed = zero(surrs[1])
    for i in eachindex(surrs[1])
       surr75[i] = quantile([surr[i] for surr in surrs], .75)
       surr05[i] = quantile([surr[i] for surr in surrs], .05)
       surr25[i] = quantile([surr[i] for surr in surrs], .25)
       surr95[i] = quantile([surr[i] for surr in surrs], .95)
       surrmed[i] = median([surr[i] for surr in surrs])
    end
    return surr05, surr25, surrmed, surr75,surr95
end

function plot_surr(dates, org, surrs)
    surr05, surr25,surrmed, surr75, surr95 = surrmetrics(org, surrs)
    #p = plot(dates, surr25, fillrange=surr75, color=:orange, line=nothing, label="Uncertainty box")
    p = plot(dates, surrs, color=:grey, alpha=0.3, label="")
    plot!(p, dates, org, color=:blue, label="")#, label="Original Imf" )
    plot!(p, dates, surrmed, color=:red, label="")#, label="Median Imf")
    plot!(p, dates, surr05, label="")#, label="5th Percentile IMF")
    plot!(p, dates, surr95,label="")#, label="95th Percentile IMF")
    p
    #return surr25, surr75
end

function plot_ceemd_uncert(dates, ts, xs)
    imfplots = []
    for i in 1:5
        org, surr = bootstrap(ts, iaaft, y->ceemdan(y, i+1)[:,end-1],200)
        p = plot_surr(dates, org, surr)
        push!(imfplots,p)
    end
    p_org = plot(dates, ts, label="")
    fullplot = plot(p_org, imfplots..., layout=(3,2))
    fullplot
end

#=
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
=#
