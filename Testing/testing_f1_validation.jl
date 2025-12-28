using IntervalArithmetic
using Roots
# using Statistics


include("functions_helpers.jl")
include("functions_stability.jl")
include("functions_validation.jl")
include("functions_f1.jl")

κ = 0.8
μ = 0.005

(θmin, θmax) = Intθ_1(κ, μ)

br = brackets_signchanges(θ -> f₁(θ, κ, μ), θmin, θmax)

root = find_zero(θ -> f₁(θ, κ, μ), br[1]; atol=1e-12)

Iκ = interval(κ, 1e-5; format=:midpoint)
Iμ = interval(μ, 1e-5; format=:midpoint)

ok, data = validate_root_sign_mono(
    f₁,
    f₁θ,
    root,
    interval(0, sup(θ1_max(Iκ, Iμ))),
    Iκ,
    Iμ;
    ntiles=16,
    min_width=1e-12,
    min_radius=1e-12,
    verbose=true,
)

complements_nonzero(
    f₁,
    data.Iθ,
    data.J,
    Iκ,
    Iμ;
    ntiles=16,
)

classify_equilibrium(data.Iexist, Iκ, Iμ)