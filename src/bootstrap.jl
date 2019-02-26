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
    surr25 = zero(org)
    surr75 = zero(org)
    surr05 = zero(org)
    surr95 = zero(org)
    surrmed = zero(org)
    for i in eachindex(org)
       surr75[i] = quantile([surr[i] for surr in surrs], .75)
       surr05[i] = quantile([surr[i] for surr in surrs], .05)
       surr25[i] = quantile([surr[i] for surr in surrs], .25)
       surr95[i] = quantile([surr[i] for surr in surrs], .95)
       surrmed[i] = median([surr[i] for surr in surrs])
    end
    #p = plot(dates, surr25, fillrange=surr75, color=:orange, line=nothing, label="Uncertainty box")
    p = plot(dates, surrs, color=:grey, alpha=0.3, label="")
    plot!(p, dates, org, color=:blue, label="Original Imf" )
    plot!(p, dates, surrmed, color=:red, label="Median Imf")
    plot!(p, dates, surr05, label="5th Percentile IMF")
    plot!(p, dates, surr95, label="95th Percentile IMF")
    p
    #return surr25, surr75
end

function plot_ceemd_uncert(dates, ts, xs)
    imfplots = []
    for i in 1:7
        org, surr = bootstrap(ts, randomshuffle, y->ceemd(y, xs,num_imfs=i)[end-1],200)
        p = plot_surr(dates, org, surr)
        push!(imfplots,p)
    end
    p_org = plot(dates, ts)
    fullplot = plot(p_org, imfplots..., layout=(4,2))
    fullplot
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
