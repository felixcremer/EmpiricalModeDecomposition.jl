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
colominas2014_s()

First testdata from Colominas 2014 et al.
http://dx.doi.org/10.1016/j.bspc.2014.06.009

"""
function colominas2014_s()
  n=1:1000
  s1 = sinpi.(2 .* 0.255 .* (n .- 501))
  s1[1:500] .= 0
  s1[751:1000] .=0
  s2 = sinpi.(2 .* 0.065 .*(n .- 1))
  return n, s1+s2, [s1,s2]
end

"""
colominas2014_x()

Second testdata from Colominas 2014 et al.
http://dx.doi.org/10.1016/j.bspc.2014.06.009
It is not exactly the same, because ϕ is not feasable because the arccos is only defined for values between -1 and 1.

"""
function colominas2014_x()
  n=1:1000
  fmax=3//32
  fmin=9//128
  ϕ = acos((fmax -fmin)/(3fmin+fmax))
  x1=3 .* exp.((.-((x .-500)./100).^2).*π) .* cospi.(5//8 .*(n.-1000))
  x2 = @. cospi((fmax+fmin)*(n - 1000)+(fmax-fmin)*500(sinpi(n/500)+ϕ-sin(ϕ)))
  x3 = exp.((.-((n.-1000)/200).^2) .*π) .* cospi.((7 .//128).*(n.-1000))
  n, x1+x2+x3, [x1,x2,x3]
end
"""
fosso2014()

Make testdata from Fosso 2014 et al.
"""
function fosso(x=0:.001:2)
  x1 = @. 0.7 * sinpi(16 * x)
  x2 = @. 0.7 * sinpi(48 * x)
  x3 = @. 1.4 * sinpi(60 * x)
  return x, x1 .+ x2 .+ x3, [x1, x2, x3]
end
