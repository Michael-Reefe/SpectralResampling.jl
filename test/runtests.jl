using SpectralResampling
using Test

@testset "SpectralResampling.jl" begin

    ### -- just testing for consistency --

    λ = 5:0.1:12
    newλ = get_logarithmic_λ(λ)[1:end-1]
    newλ2 = get_logarithmic_λ(λ, log(newλ[2]/newλ[1]))
    @test all(newλ .≈ newλ2)
    @test all(get_linear_λ(get_logarithmic_λ(λ)) .≈ λ)

    # 1D test
    F = randn(length(λ))
    e = 0.1 .* F
    m = falses(size(F))
    m[5] = true
    newF, newe, newm = resample_conserving_flux(newλ, λ, F, e, m)
    newF2, newe2 = resample_conserving_flux(newλ, λ, F, e)
    newF3 = resample_conserving_flux(newλ, λ, F)
    @test all(newF .== newF2 .== newF3)
    @test all(newe .== newe2)

    # 2D test
    F = randn(length(λ), 2)
    e = 0.1 .* F
    m = falses(size(F))
    m[:,2] .= true
    newF, newe, newm = resample_conserving_flux(newλ, λ, F, e, m)
    newF2, newe2 = resample_conserving_flux(newλ, λ, F, e)
    newF3 = resample_conserving_flux(newλ, λ, F)
    @test all(newF .== newF2 .== newF3)
    @test all(newe .== newe2)

    ## -- testing for correctness --

    λ = [5., 5.1, 5.2, 5.3, 5.4]
    newλ = [5.16, 5.35]
    F = [1., 1.2, 1.4, 1.2, 1.]
    e = 0.1 .* F
    m = BitVector([0, 0, 1, 0, 0])
    newF, newe, newm = resample_conserving_flux(newλ, λ, F, e, m)
    @test all(newF .≈ [1.3052631578947367, 1.1000000000000003])
    @test all(newe .≈ [0.09122132228755071, 0.0781024967590666])
    @test all(newm .== BitVector([1, 0]))

end
