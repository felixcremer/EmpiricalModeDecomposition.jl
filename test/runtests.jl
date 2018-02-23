using EmpiricalModeDecomposition
using Base.Test

# write your own tests here
@testset "ismonotonic" begin
    @test EmpiricalModeDecomposition.ismonotonic(1:100)
    @test EmpiricalModeDecomposition.ismonotonic(100:-1:1)
    @test EmpiricalModeDecomposition.ismonotonic([1.0,2.0,3.0,6.0])
    @test EmpiricalModeDecomposition.ismonotonic([6.0,5.0,1.0])
    @test EmpiricalModeDecomposition.ismonotonic([1.0,3.0,2.0]) == false
    @test EmpiricalModeDecomposition.ismonotonic([3.0,1.0,2.0]) == false
end

@testset "emd" begin
    x = -1:0.1:2π+1
    measurements = sin.(x)
    imfs = emd(measurements,x)
    @test isapprox(imfs[1], measurements, rtol=0.001)
end

@testset "localmaxmin" begin
    maxes = Int[]
    mins  = Int[]
    x = zeros(100)
    x_max = 2:10:100
    x[x_max] = 1
    x_min = 5:10:100
    x[x_min]=-1
    EmpiricalModeDecomposition.localmaxmin!(x,maxes,mins)
    @test maxes == collect(x_max)
    @test mins  == collect(x_min)
end
