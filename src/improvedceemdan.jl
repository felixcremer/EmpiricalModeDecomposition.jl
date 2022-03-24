"""
iceemdan(measurements, times)

    Decompose a time series of measurements with time steps times using the improved complete ensemble empirical mode decomposition with
    adaptive noise (ICEEMDAN). 
    The results are a vector of intrinsic mode decompositions. 
    This follows:
    [1]M. A. Colominas, G. Schlotthauer, und M. E. Torres, „Improved complete ensemble EMD: A suitable tool for biomedical signal processing“, Biomedical Signal Processing and Control, Bd. 14, S. 19–29, Nov. 2014, doi: 10.1016/j.bspc.2014.06.009.
"""
function iceemdan(measurements, times; num_imfs=log2(length(measurements)), numtrails=100)
    eps0 = 0.2
    β_0 = eps0 * std(measurements)
    whitenoise_ens = [randn(length(measurements)) for i in 1:numtrails]
    whitenoise_imf_ens = emd.(whitenoise_ens, Ref(times))
    imfnum = min(minimum(length.(whitenoise_imf_ens)), num_imfs)
    measurements_ens = [measurements + β_0*wnimf[1] for wnimf in whitenoise_imf_ens]
    prevresidual = mean((measurements,) .- sift.(measurements_ens))
    nextresidual = similar(prevresidual)
    imfs = Vector{eltype(measurements)}[measurements - prevresidual]
    #@show typeof(imfs)
    for i in 1:imfnum
        beta = eps0 * std(prevresidual) 
        #@show beta
        measurements_ens .= [prevresidual .+ beta .* wnimf[i] for wnimf in whitenoise_imf_ens]
        nextresidual .= mean((prevresidual,) .- sift.(measurements_ens))
        push!(imfs, prevresidual .- nextresidual)
        prevresidual .= nextresidual
    end
    push!(imfs, prevresidual)
    imfs

end