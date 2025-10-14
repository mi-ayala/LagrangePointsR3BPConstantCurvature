using IntervalArithmetic
using Plots
using PlotlyJS
using Roots
plotlyjs()

include("functions.jl")

κ = 0.5
μ = 0.02


interior_linspace(lo, hi, n; margin=0.1) = begin
    a = lo + margin * (hi - lo)
    b = hi - margin * (hi - lo)
    range(a, b; length=n)
end

theta_range = interior_linspace(Intθ_1(κ, μ)[1], Intθ_1(κ, μ)[2], 100)

# f₁(θ, κ, μ)
# f₁θ(θ, κ, μ)

Plots.plot(theta_range, @. f₁(theta_range, κ, μ))


### A single root example

### Estimate by sign changes
for i in range(Intθ_1(κ, μ)[1], Intθ_1(κ, μ)[2]; length=400)
    println("f₁($i) = ", f₁(i, κ, μ))
end
br = brackets_f₁(κ, μ; n=800)
root = find_zero(θ -> f₁(θ, κ, μ), br[1])

delta = 0.0001

f₁(interval(root + delta), interval(κ), interval(μ))
f₁(interval(root - delta), interval(κ), interval(μ))

f₁θ(interval(root, delta; format=:midpoint), interval(κ), interval(μ))


parameter_delta = 0.00001

f₁(interval(root + delta), interval(κ, parameter_delta; format=:midpoint), interval(μ, parameter_delta; format=:midpoint))

f₁(interval(root - delta), interval(κ, parameter_delta; format=:midpoint), interval(μ, parameter_delta; format=:midpoint))

f₁θ(interval(root, delta; format=:midpoint), interval(κ, parameter_delta; format=:midpoint), interval(μ, parameter_delta; format=:midpoint))


function θ_max(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    y = tan(r * q₂(κ, μ))     # numerator
    x = one(μ) + μ            # denominator (kept as same type)
    return atan(y, x)
end



f₁θ(interval(0, θ_max(interval(κ), interval(μ))), interval(κ), interval(μ))


# # bounds for the box
# μT = 0.0886765
# κ_lo, κ_hi = 0.0, (0.5 * pi)^2
# μ_lo, μ_hi = 0.0, μT


