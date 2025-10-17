using IntervalArithmetic
using Plots
using PlotlyJS
using Roots
plotlyjs()

include("functions.jl")

κ = 0.5
μ = 0.02

interior_linspace(lo, hi, n; margin=0.01) = begin
    a = lo + margin * (hi - lo)
    b = hi - margin * (hi - lo)
    range(a, b; length=n)
end

theta_range = interior_linspace(Intθ_1(κ, μ)[1], Intθ_1(κ, μ)[2], 100)

Plots.plot(
    Plots.plot(theta_range, @. f₁(theta_range, κ, μ)),
    Plots.plot(theta_range, @. f₁θ(theta_range, κ, μ)))


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


f₁θ(interval(0, θ1_max(interval(κ), interval(μ))), interval(κ), interval(μ))


# # bounds for the box
# μT = 0.0886765
# κ_lo, κ_hi = 0.0, (0.5 * pi)^2
# μ_lo, μ_hi = 0.0, μT


# example
κ = 0.5
μ = 0.02
r_κ = 0.001
r_μ = 0.001

validate_f1_root(root, κ, μ, r_κ, r_μ; verbose=true)

function push_kappa_interval(κ_center, μ_center, root;
    Δκ=1e-3, κmax=(π / 2)^2,
    r_κ=1e-8, r_μ=1e-8,
    verbose=false)

    # first validation
    ok = validate_f1_root(root, κ_center, μ_center, r_κ, r_μ; verbose)
    ok || error("initial point not valid")

    κ_best, root_best = κ_center, root

    while κ_center < κmax
        # next center value (float)
        κ_center_next = κ_center + Δκ

        # build intervals for rigorous check
        Iκ = interval(κ_center_next, r_κ; format=:midpoint)
        Iμ = interval(μ_center, r_μ; format=:midpoint)

        # compute the root at midpoint values (for seeding)
        root_next = find_zero(θ -> f₁(θ, mid(Iκ), mid(Iμ)), root)

        if !isfinite(root_next)
            verbose && println("nonfinite root at κ=$(mid(Iκ))")
            break
        end

        if validate_f1_root(root_next, mid(Iκ), mid(Iμ), r_κ, r_μ; verbose)
            κ_center, root = mid(Iκ), root_next
            κ_best, root_best = κ_center, root
            verbose && println("✓ valid at κ=$(round(κ_center, sigdigits=6))")
        else
            verbose && println("✗ failed at κ=$(round(mid(Iκ), sigdigits=6))")
            break
        end
    end

    return interval(κ_best, r_κ; format=:midpoint), root_best
end


κ0, μ = 0.5, 0.02

Iκ_final, root_final = push_kappa_interval(κ0, μ, root; Δκ=1e-3, r_κ=1e-8, verbose=true)


θc = root
rθ = interval(0, θ1_max(interval(κ), interval(μ)))
κ = 0.5
μ = 0.02
rκ = 1e-8
rμ = 1e-8

status, Iθ_cert, dI_cert = certify_deriv_sign_centered(θc, κ, μ, rκ, rμ; verbose=true)


Iκ = interval(κ, rκ; format=:midpoint)
Iμ = interval(μ, rμ; format=:midpoint)

Iθ = interval(0, θ1_max(Iκ, Iμ))

Jθ = Iθ_cert
κ, μ = 0.5, 0.02
rκ, rμ = 1e-8, 1e-8

ok = validate_complements(Iθ, Jθ, κ, μ, rκ, rμ; verbose=true)
