"""
    maketestdata(seed)

Return a simple example time series with the composing parts.
"""
function maketestdata(seed=1234, nyears=10,step=12)
  Random.seed!(seed)
  ## simulate data of length....
  NpY   = 365/step         # samples/year
  t = 1:step:(nyears*365)
  #t     = t./NpY # your time vector

  # constant seasonal cycle
  A   = 2          # amplitude
  phi = 13 * pi ./12   # initial phase
  S   = A .* cos.(2 .*pi .* t ./ 365   .+ phi)


  # generate a linear trend
  T = 0.1 .+ 0.02.*t
  @debug T

  # some other oscillation
  a = 0.2
  b = 0.1
  C = 1 .+ a .*cos.(b*2*pi*t./NpY)

  # simple (coloured) noise
  φ = 0.3 # strengh of autocorrelation in noise
  E = randn(length(t)) .* 0.1

  for i = 2:length(t)
    E[i] = φ*E[i-1]+(1-φ)*E[i]
  end

  X = @. S * C + 2 * E + T
  t, X, [S, C, E, T]
end

"""
colominas2014_s()

First test data from Colominas 2014 et al.
http://dx.doi.org/10.1016/j.bspc.2014.06.009

"""
function colominas2014_s()
  n = 1:1000
  s1 = sinpi.(2 .* 0.255 .* (n .- 501))
  s1[1:500] .= 0
  s1[751:1000] .=0
  s2 = sinpi.(2 .* 0.065 .*(n .- 1))
  n, s1+s2, [s1,s2]
end

"""
    colominas2014_x()

Second test data from Colominas 2014 et al.
http://dx.doi.org/10.1016/j.bspc.2014.06.009
It is not exactly the same, because ϕ is not feasable since the arccos is only
defined for values between -1 and 1.
"""
function colominas2014_x()
  n=1:2000
  fmax=3//32
  fmin=9//128
  ϕ = -acos((fmax -fmin)/(3fmin+fmax))
  x1=3 .* exp.((.-((n .-500)./100).^2).*π) .* cospi.(5//8 .*(n.-1000))
  x2 = cos.(π .* (fmax .+ fmin) .* (n .- 1000) .+ (fmax .- fmin) .* 500 .* (sinpi.(n ./ 500) .+ ϕ .- sin.(ϕ)))
  x3 = exp.((.-((n.-1000)/200).^2) .*π) .* cospi.((7 .//128).*(n.-1000))
  n, x1+x2+x3, [x1,x2,x3]
end

"""
    fosso2014()

Make testdata from Fosso 2017 et al.
"""
function fosso(x=0:.001:2)
  x1 = @. 0.7 * sinpi(16 * x)
  x2 = @. 0.7 * sinpi(48 * x)
  x3 = @. 1.4 * sinpi(60 * x)
  x, x1 .+ x2 .+ x3, [x1, x2, x3]
end

function sawtooth()
  a1 = 35.:-1:6.
  a3 = 33.:-1:4
  a2 = 27.:-1:-2
  @show length.([a1,a2,a3])
  b = 2.:-1:-2.
  @show length(b)
  c = 29:-1:25.
  @show length(c)
  sawtooth = append!(collect(b), a1,a2,a3,c) 
  msaw = mean(sawtooth)
  return 1:100, sawtooth, [sawtooth .- msaw, zero(sawtooth) .+ msaw]
end