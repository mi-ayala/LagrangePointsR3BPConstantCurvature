using IntervalArithmetic
using Plots
using PlotlyJS
using Roots
using Statistics
plotlyjs()

include("functions.jl")
include("validation_functions.jl")

κ = 0.1
μ = 0.01

interior_linspace(lo, hi, n; margin=0.01) = begin
    a = lo + margin * (hi - lo)
    b = hi - margin * (hi - lo)
    range(a, b; length=n)
end

theta_range = interior_linspace(Intθ_4(κ, μ)[1], Intθ_4(κ, μ)[2], 100)

### Plot function and derivative
Plots.plot(
    Plots.plot(theta_range, @. f₄(theta_range, κ, μ)),
    Plots.plot(theta_range, @. f₄θ(theta_range, κ, μ)))

### A single root example
### Estimate by sign changes
for i in range(Intθ_4(κ, μ)[1], Intθ_4(κ, μ)[2]; length=400)
    println("f₄($i) = ", f₄(i, κ, μ))
end
br = brackets_f₄(κ, μ; n=800)
root = find_zero(θ -> f₄(θ, κ, μ), br[1])

delta = 0.00001


f₄(interval(root + delta), interval(κ), interval(μ))
f₄(interval(root - delta), interval(κ), interval(μ))
f₄θ(interval(root, delta; format=:midpoint), interval(κ), interval(μ))

parameter_delta = 0.0000001
f₄(interval(root + delta), interval(κ, parameter_delta; format=:midpoint), interval(μ, parameter_delta; format=:midpoint))

f₄(interval(root - delta), interval(κ, parameter_delta; format=:midpoint), interval(μ, parameter_delta; format=:midpoint))

f₄θ(interval(root, delta; format=:midpoint), interval(κ, parameter_delta; format=:midpoint), interval(μ, parameter_delta; format=:midpoint))

f₄θ(interval(inf(θ4_min(interval(κ), interval(μ))), sup(θ4_max(interval(κ), interval(μ)))), interval(κ), interval(μ))

### Testing validation function
parameter_delta_kappa = 0.00000001
parameter_delta_mu = 0.00000001
validate_f4_root(root, κ, μ, parameter_delta_kappa,
    parameter_delta_mu; verbose=true)


# refine_root_interval(root, κ, μ, parameter_delta_kappa,
#     parameter_delta_mu,interval(0.226063, 0.31039) )

# validate_and_classify(root, κ, μ,
#     parameter_delta_kappa,
#     parameter_delta_mu; verbose=true)

# Iκ = interval(κ, parameter_delta_kappa; format=:midpoint)
# Iμ = interval(μ, parameter_delta_mu; format=:midpoint)
# classify_equilibrium(interval(0.268185, 0.268268), Iκ, Iμ; verbose=true)