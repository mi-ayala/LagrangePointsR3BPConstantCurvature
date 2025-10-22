using IntervalArithmetic
using Plots
using PlotlyJS
using Roots
plotlyjs()

include("functions.jl")
include("validation_functions.jl")

κ = 2.0
μ = 0.05

interior_linspace(lo, hi, n; margin=0.01) = begin
    a = lo + margin * (hi - lo)
    b = hi - margin * (hi - lo)
    range(a, b; length=n)
end

theta_range = interior_linspace(Intθ_1(κ, μ)[1], Intθ_1(κ, μ)[2], 100)

### Plot function and derivative
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

delta = 0.00001

f₁(interval(0.901893), interval(κ), interval(μ))

f₁(interval(root + delta), interval(κ), interval(μ))
f₁(interval(root - delta), interval(κ), interval(μ))
f₁θ(interval(root, delta; format=:midpoint), interval(κ), interval(μ))

parameter_delta = 0.0000001
f₁(interval(root + delta), interval(κ, parameter_delta; format=:midpoint), interval(μ, parameter_delta; format=:midpoint))

f₁(interval(root - delta), interval(κ, parameter_delta; format=:midpoint), interval(μ, parameter_delta; format=:midpoint))

f₁θ(interval(root, delta; format=:midpoint), interval(κ, parameter_delta; format=:midpoint), interval(μ, parameter_delta; format=:midpoint))

f₁θ(interval(0, θ1_max(interval(κ), interval(μ))), interval(κ), interval(μ))

### Testing validation function
parameter_delta_kappa = 0.0000001
parameter_delta_mu = 0.0000001
validate_f1_root(root, κ, μ, parameter_delta_kappa,
    parameter_delta_mu; verbose=true)

# Maximize rμ so that validate_f1_root(root, κ, μ, rκ_fixed, rμ) == true
function maximize_delta_mu(root, κ, μ,
    rκ_fixed;
    rμ0=1e-8,       # starting radius for μ
    rμ_max=1.0,     # hard cap
    grow=1.8,       # multiplicative growth
    tol=1e-12,      # final bisection tolerance
    max_steps=200,
    verbose=false)

    # helper: predicate we’re maximizing against
    ok(rμ) = validate_f1_root(root, κ, μ, rκ_fixed, rμ; verbose=false)

    # normalize types
    T = promote_type(typeof(root), typeof(κ), typeof(μ), typeof(rκ_fixed),
        typeof(rμ0), typeof(rμ_max), typeof(tol))
    rμ = T(rμ0)
    rμ_cap = T(rμ_max)
    tolT = T(tol)

    # 0) ensure we start from a valid radius (shrink if needed)
    steps = 0
    while steps < max_steps && !ok(rμ) && rμ > tolT
        rμ /= 2
        steps += 1
    end
    if !ok(rμ)
        verbose && println("No valid μ-radius found (even below tol).")
        return (false, κ, μ, rκ_fixed, zero(T))
    end

    # 1) grow until it fails or hits rμ_cap
    last_good = rμ
    first_bad = T(NaN)
    while steps < max_steps && rμ < rμ_cap
        rμ_try = min(rμ * T(grow), rμ_cap)
        if ok(rμ_try)
            last_good = rμ = rμ_try
            verbose && println("✓ rμ = ", rμ)
        else
            first_bad = rμ_try
            verbose && println("✗ fail at rμ = ", rμ_try)
            break
        end
        steps += 1
    end

    # 2) if we never failed, we’re done (hit cap)
    if !isfinite(first_bad)
        return (true, κ, μ, rκ_fixed, last_good)
    end

    # 3) bisect between last_good and first_bad
    lo, hi = last_good, first_bad
    while (hi - lo) > tolT && steps < max_steps
        mid = (lo + hi) / 2
        if ok(mid)
            lo = mid
        else
            hi = mid
        end
        steps += 1
    end
    return (true, κ, μ, rκ_fixed, lo)
end


function maximize_delta_kappa(root, κ, μ,
    rμ_fixed;
    rκ0::Real=1e-8,     # starting radius for κ
    rκ_max::Real=1.0,   # hard cap you allow
    grow::Real=1.8,     # multiplicative growth
    tol::Real=1e-12,   # final bisection tolerance
    max_steps::Int=200,
    verbose::Bool=false)

    # predicate: does the validator still certify with this κ-radius?
    ok(rκ) = validate_f1_root(root, κ, μ, rκ, rμ_fixed; verbose=false)

    # normalize numeric type
    T = promote_type(typeof(root), typeof(κ), typeof(μ),
        typeof(rμ_fixed), typeof(rκ0), typeof(rκ_max), typeof(tol))
    rκ = T(rκ0)
    rκ_cap = T(rκ_max)
    tolT = T(tol)
    gT = T(grow)

    # 0) ensure starting point is valid (shrink if needed)
    steps = 0
    while steps < max_steps && !ok(rκ) && rκ > tolT
        rκ /= 2
        steps += 1
    end
    if !ok(rκ)
        verbose && println("No valid κ-radius found (even below tol).")
        return (false, κ, μ, zero(T), rμ_fixed)
    end

    # 1) grow until first failure or cap
    last_good = rκ
    first_bad = T(NaN)
    while steps < max_steps && rκ < rκ_cap
        rκ_try = min(rκ * gT, rκ_cap)
        if ok(rκ_try)
            last_good = rκ = rκ_try
            verbose && println("✓ rκ = ", rκ)
        else
            first_bad = rκ_try
            verbose && println("✗ fail at rκ = ", rκ_try)
            break
        end
        steps += 1
    end

    # 2) if never failed, we hit the cap
    if !isfinite(first_bad)
        return (true, κ, μ, last_good, rμ_fixed)
    end

    # 3) bisect between last_good and first_bad
    lo, hi = last_good, first_bad
    while (hi - lo) > tolT && steps < max_steps
        mid = (lo + hi) / 2
        if ok(mid)
            lo = mid
        else
            hi = mid
        end
        steps += 1
    end

    return (true, κ, μ, lo, rμ_fixed)
end



### Testing validation function
parameter_delta_kappa = 0.000001
parameter_delta_mu = 0.000001
validate_f1_root(root, κ, μ, parameter_delta_kappa,
    parameter_delta_mu; verbose=true)


ok, κc, μc, rκ_star, rμ = maximize_delta_mu(root, κ, μ,
    parameter_delta_kappa;
    rμ0=parameter_delta_mu,   # starting radius for μ
    rμ_max=0.0886765,                  # hard cap (set as you like)
    grow=1.8,                   # multiplicative growth
    tol=1e-12,             # final bisection tolerance
    max_steps=200,
    verbose=true)


validate_f1_root(root, κc, μc, rκ_star, rμ; verbose=true)

ok, κc, μc, rκ, rμ = maximize_delta_kappa(root, κ, μ, rμ;
    rκ0=rκ_star,     # starting radius for κ
    rκ_max=(0.5 * pi)^2,   # hard cap 
    grow=1.8,     # multiplicative growth
    tol=1e-12,   # final bisection tolerance
    max_steps=200,
    verbose=true)


"""
    maximize_box(root, κc, μc; rκ0=1e-8, rμ0=1e-8,
                 κ_hi=(0.5*π)^2, μT=0.0886765,
                 grow=1.8, tol=1e-12, itmax=20, verbose=false)

Alternating maximization of (rκ, rμ) around (κc, μc), stopping when neither radius improves.
Returns (ok, κc, μc, rκ, rμ).
"""
function maximize_box(root, κc, μc;
    rκ0=1e-9, rμ0=1e-9,
    κ_hi=(0.5 * π)^2, μT=0.0886765,
    grow=1.2, tol=1e-12, itmax=20, verbose=false)

    rκ, rμ = rκ0, rμ0
    for t in 1:itmax
        # cap by remaining distance to boundary
        rκ_cap = min((κ_hi - κc), κ_hi)  # mild cap; tweak as you wish
        rμ_cap = min((μT - μc), μT)

        # κ step (hold rμ)
        okκ, _, _, rκ_new, _ = maximize_delta_kappa(root, κc, μc,
            rμ; rκ0=rκ, rκ_max=rκ_cap, grow=grow, tol=tol, verbose=false)

        # μ step (hold rκ)
        okμ, _, _, _, rμ_new = maximize_delta_mu(root, κc, μc,
            rκ_new; rμ0=rμ, rμ_max=rμ_cap, grow=grow, tol=tol, verbose=false)

        # stop if nothing improved meaningfully
        if (!okκ && !okμ) || (abs(rκ_new - rκ) ≤ 1e-3 * max(rκ, 1e-30) &&
                              abs(rμ_new - rμ) ≤ 1e-3 * max(rμ, 1e-30))
            rκ, rμ = max(rκ, 0.0), max(rμ, 0.0)
            ok = rκ > 0 && rμ > 0
            return (ok, κc, μc, rκ, rμ)
        end

        rκ, rμ = rκ_new, rμ_new
    end
    (true, κc, μc, rκ, rμ)
end



"""
    cover_rectangle(κ_hi=(0.5*π)^2, μT=0.0886765;
        nx=12, ny=12, rκ0=1e-8, rμ0=1e-8, verbose=true)

Coarse seeds + local box maximization, returns a vector of (ok, κc, μc, rκ, rμ)
suitable for plot_validated_boxes().
"""
function cover_rectangle(κ_hi=(0.5 * π)^2, μT=0.0886765;
    nx=200, ny=200, rκ0=1e-9, rμ0=1e-9, verbose=true)

    κs = range(0.02, κ_hi; length=nx)
    μs = range(0.0, μT; length=ny)

    boxes = Tuple{Bool,Float64,Float64,Float64,Float64}[]

    for κc in κs, μc in μs
        # find a root to seed the local box (skip if bracketing fails)
        br = brackets_f₁(κc, μc; n=200)
        if isempty(br)
            # try a neutral theta guess if you like; otherwise skip
            verbose && @debug "No bracket at (κ,μ)=($κc,$μc); skipping seed."
            continue
        end
        root = find_zero(θ -> f₁(θ, κc, μc), br[1])

        ok, κc2, μc2, rκ, rμ = maximize_box(root, κc, μc;
            rκ0=rκ0, rμ0=rμ0, κ_hi=κ_hi, μT=μT, verbose=false)
        ok && push!(boxes, (ok, κc2, μc2, rκ, rμ))
    end
    boxes
end


boxes = cover_rectangle()
plot_validated_boxes(boxes; color=:teal, alpha=0.45, title="Validated cover of (κ, μ)")
