using Test
using ForwardDiff
using IntervalArithmetic

include("../src/functions_helpers.jl")
include("../src/functions_f1.jl")
include("../src/functions_f2.jl")
include("../src/functions_f3.jl")
include("../src/functions_f4.jl")
include("../src/functions_g.jl")
include("../src/functions_gneg.jl")
include("../src/functions_f1neg.jl")
include("../src/functions_f2neg.jl")
include("../src/functions_f3neg.jl")
include("../src/functions_f4neg.jl")

const ATOL = 1e-10
const RTOL = 1e-8

function check_derivative(label, ad, an; atol=ATOL, rtol=RTOL)
    println()
    println("=== $label ===")
    @show an ad
    diff = an - ad
    @show diff
    @test isapprox(an, ad; atol=atol, rtol=rtol)
end

@testset "Analytic derivatives vs ForwardDiff" begin
    θ0 = 0.35
    x0 = 0.35
    κ0 = 0.7
    κneg0 = -0.7
    μ0 = 0.02

    @testset "Positive κ functions" begin
        check_derivative(
            "f₁θ",
            ForwardDiff.derivative(t -> f₁(t, κ0, μ0), θ0),
            f₁θ(θ0, κ0, μ0),
        )

        check_derivative(
            "f₂θ",
            ForwardDiff.derivative(t -> f₂(t, κ0, μ0), θ0),
            f₂θ(θ0, κ0, μ0),
        )

        check_derivative(
            "f₂θθ",
            ForwardDiff.derivative(t -> f₂θ(t, κ0, μ0), θ0),
            f₂θθ(θ0, κ0, μ0),
        )

        check_derivative(
            "f₂κ",
            ForwardDiff.derivative(k -> f₂(θ0, k, μ0), κ0),
            f₂κ(θ0, κ0, μ0),
        )

        check_derivative(
            "f₂μ",
            ForwardDiff.derivative(m -> f₂(θ0, κ0, m), μ0),
            f₂μ(θ0, κ0, μ0),
        )

        check_derivative(
            "f₃θ",
            ForwardDiff.derivative(t -> f₃(t, κ0, μ0), θ0),
            f₃θ(θ0, κ0, μ0),
        )

        check_derivative(
            "f₃θθ",
            ForwardDiff.derivative(t -> f₃θ(t, κ0, μ0), θ0),
            f₃θθ(θ0, κ0, μ0),
        )

        check_derivative(
            "f₄θ",
            ForwardDiff.derivative(t -> f₄(t, κ0, μ0), θ0),
            f₄θ(θ0, κ0, μ0),
        )

        check_derivative(
            "gθ",
            ForwardDiff.derivative(t -> g(t, κ0, μ0), x0),
            gθ(x0, κ0, μ0),
        )
    end

    @testset "Negative κ functions" begin
        check_derivative(
            "f₁negθ",
            ForwardDiff.derivative(t -> f₁neg(t, κneg0, μ0), θ0),
            f₁negθ(θ0, κneg0, μ0),
        )

        check_derivative(
            "f₂negθ",
            ForwardDiff.derivative(t -> f₂neg(t, κneg0, μ0), θ0),
            f₂negθ(θ0, κneg0, μ0),
        )

        check_derivative(
            "f₂negθθ",
            ForwardDiff.derivative(t -> f₂negθ(t, κneg0, μ0), θ0),
            f₂negθθ(θ0, κneg0, μ0),
        )

        check_derivative(
            "f₂negκ",
            ForwardDiff.derivative(k -> f₂neg(θ0, k, μ0), κneg0),
            f₂negκ(θ0, κneg0, μ0),
        )

        check_derivative(
            "f₂negμ",
            ForwardDiff.derivative(m -> f₂neg(θ0, κneg0, m), μ0),
            f₂negμ(θ0, κneg0, μ0),
        )

        check_derivative(
            "f₃negθ",
            ForwardDiff.derivative(t -> f₃neg(t, κneg0, μ0), θ0),
            f₃negθ(θ0, κneg0, μ0),
        )

        check_derivative(
            "f₃negθθ",
            ForwardDiff.derivative(t -> f₃negθ(t, κneg0, μ0), θ0),
            f₃negθθ(θ0, κneg0, μ0),
        )

        check_derivative(
            "f₄negθ",
            ForwardDiff.derivative(t -> f₄neg(t, κneg0, μ0), θ0),
            f₄negθ(θ0, κneg0, μ0),
        )

    end
end
