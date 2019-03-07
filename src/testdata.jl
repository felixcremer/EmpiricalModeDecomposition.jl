"""
maketestdata(seed)

Return a simple example time series with the composing parts
"""
function maketestdata(seed)
  Random.seed!(seed)
  ## simulate data of length....
  N_tim = 240
  NpY   = 24         # samples/year
  t     = 1:N_tim
  t     = t./NpY # your time vector

  # constant seasonal cycle
  A   = 2          # amplitude
  phi = 13 * pi/12   # initial phase
  S   = A .* cos.(2 .* pi .* t .+ phi)

  # generate a linear trend
  T = 0.1 .+ 0.2.*t
  @show T

  # some other oscillation
  a = 0.2
  b = 0.1
  C = 1 .+ a .*cos.(b*2*pi*t)

  # simple (coloured) noise
  φ = 0.3 # strengh of autocorrelation in noise
  E = randn(N_tim) .* 0.1

  for i = 2:N_tim
    E[i] = φ*E[i-1]+(1-φ)*E[i]
  end

  X = S .* C .+ 2. .*E .+ T
  t, X, [S,C,E,T]
end

"""
colominas()

Make testdata from Colominas 2014 et al.
"""
function colominas()
  x = 1:1000
  s1 = sinpi.(2 .* 0.255 .* (x .- 501))
  s1[1:500] .= 0
  s1[751:1000] .=0
  s2 = sinpi.(2 .* 0.065 .*(x .- 1))
  return x, s1+s2, (s1,s2)
end
