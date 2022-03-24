using EmpiricalModeDecomposition
import EmpiricalModeDecomposition: ismonotonic, localmaxmin!, get_edgepoint,
    SiftIterable
using Test
using Random

@testset "emd" begin
    x = 0:0.1:10
    measurements = sinpi.(x)
    imfs = emd(measurements,x)
    @test isapprox(imfs[1], measurements, rtol=0.001)
end
@testset "EmpiricalModeDecomposition.jl" begin
    @testset "ismonotonic" begin
        @test ismonotonic(1:100)
        @test ismonotonic(100:-1:1)
        @test ismonotonic([1.0,2.0,3.0,6.0])
        @test ismonotonic([6.0,5.0,1.0])
        @test !ismonotonic([1.0,3.0,2.0])
        @test !ismonotonic([3.0,1.0,2.0])
    end

    @testset "emd" begin
        x = -1:0.1:2π+1
        measurements = sin.(x)
        imfs = emd(measurements, x)
        @test isapprox(imfs[1], measurements, rtol=0.001)
    end

    @testset "localmaxmin" begin
        maxes = Int[]
        mins  = Int[]
        x = zeros(100)
        x_max = 2:10:100
        x[x_max] .= 1
        x_min = 5:10:100
        x[x_min] .= -1
        localmaxmin!(x, maxes, mins)
        @test maxes == collect(x_max)
        @test mins  == collect(x_min)
    end

    @testset "get_edgepoint" begin
        t = 0:0.5:10
        y = zero(t)
        y[3] = 1
        y[5] = 0.5
        y[17] = 0.5
        y[19] = 1.0
        maxes = Int[]
        mins = Int[]
        localmaxmin!(y, maxes, mins)
        @test get_edgepoint(y, t, maxes, first, !isless) == 1.5
        @test get_edgepoint(y, t, maxes, last, !isless)  == 1.5
    end

    @testset "SiftIterable" begin
        x = -1:0.1:2π+1
        measurements = sin.(x)
        imf = zero(measurements)
        for sift in Base.Iterators.take(SiftIterable(measurements,x,6),10)
            imf = sift.yvec
        end
        @test imf ≈ measurements
    end

@testset "EEMD" begin
    x = -1:0.1:2π+1
    measurements = sin.(x) .+ cos.(2*x) .+ 2 .*x
    imf = eemd(measurements, x)
    @test sum(imf) ≈ measurements
end

function maketestdata(seed)
  Random.seed!(seed)
  ## simulate data of length....
  N_tim = 240
  NpY   = 24         # samples/year
  t     = 1:N_tim
  t     = t./NpY # your time vector

  # constant seasonal cycle
  A   = 2          # amplitude
  phi = 13*pi/12   # initial phase
  S   = A*cos.(2 * pi * t+phi)

  # generate a linear trend
  T = 0.1 + 0.0002 .* t

  # some other oscillation
  a = 0.2
  b = 0.1
  C = 1+a*cos.(b*2*pi*t)

  # simple (coloured) noise
  φ = 0.3 # strengh of autocorrelation in noise
  E = randn(N_tim).*0.1

  for i = 2:N_tim
    E[i] = φ*E[i-1]+(1-φ)*E[i]
  end

  X = @. S*C + 2*E
  t, X, E
end
