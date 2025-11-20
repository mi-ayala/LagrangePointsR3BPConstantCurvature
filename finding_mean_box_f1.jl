using IntervalArithmetic
using Plots
using PlotlyJS
using Roots
using Statistics
plotlyjs()

include("functions.jl")
include("validation_functions.jl")



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

function maximize_box(root, κc, μc;
    rκ0=1e-9, rμ0=1e-9,
    κ_hi=(0.5 * π)^2, μT=0.0886765,
    grow=1.2, tol=1e-12, itmax=20,
    thicken_rounds::Int=1,     # how many "halve-big, grow-other" passes
    polish_iters::Int=2,       # small post-polish (alt κ/μ)
    verbose=false)

    # ---------- 1) baseline: your original alternating growth ----------
    rκ, rμ = rκ0, rμ0
    for _ in 1:itmax
        rκ_cap = κ_hi - κc
        rμ_cap = μT - μc

        okκ, _, _, rκ_new, _ = maximize_delta_kappa(root, κc, μc, rμ;
            rκ0=rκ, rκ_max=rκ_cap, grow=grow, tol=tol, verbose=verbose)

        okμ, _, _, _, rμ_new = maximize_delta_mu(root, κc, μc, rκ_new;
            rμ0=rμ, rμ_max=rμ_cap, grow=grow, tol=tol, verbose=verbose)

        # stop if no meaningful improvement
        if (!okκ && !okμ) ||
           (abs(rκ_new - rκ) ≤ 1e-3 * max(rκ, 1e-30) &&
            abs(rμ_new - rμ) ≤ 1e-3 * max(rμ, 1e-30))
            rκ, rμ = max(rκ, 0.0), max(rμ, 0.0)
            break
        end
        rκ, rμ = rκ_new, rμ_new
    end

    # If either radius stayed zero, we can't thicken
    if !(rκ > 0 && rμ > 0)
        return (false, κc, μc, rκ, rμ)
    end

    # ---------- 2) thickening: halve the bigger, grow the smaller ----------
    for pass in 1:thicken_rounds
        rκ_cap = κ_hi - κc
        rμ_cap = μT - μc

        if rκ ≥ rμ
            rκ_half = max(rκ * 0.5, 1e-30)
            okμ, _, _, _, rμ_new = maximize_delta_mu(root, κc, μc, rκ_half;
                rμ0=rμ, rμ_max=rμ_cap, grow=grow, tol=tol, verbose=verbose)
            if okμ && rμ_new > rμ
                rκ, rμ = rκ_half, rμ_new
            else
                # even if μ didn't grow, keep the halved κ only if it helps thickness
                rκ = rκ_half
            end
        else
            rμ_half = max(rμ * 0.5, 1e-30)
            okκ, _, _, rκ_new, _ = maximize_delta_kappa(root, κc, μc, rμ_half;
                rκ0=rκ, rκ_max=rκ_cap, grow=grow, tol=tol, verbose=verbose)
            if okκ && rκ_new > rκ
                rκ, rμ = rκ_new, rμ_half
            else
                rμ = rμ_half
            end
        end
    end

    # ---------- 3) tiny polish: a couple of alternations ----------
    for _ in 1:polish_iters
        rκ_cap = κ_hi - κc
        rμ_cap = μT - μc

        okκ, _, _, rκ_new, _ = maximize_delta_kappa(root, κc, μc, rμ;
            rκ0=rκ, rκ_max=rκ_cap, grow=grow, tol=tol, verbose=verbose)
        okμ, _, _, _, rμ_new = maximize_delta_mu(root, κc, μc, rκ_new;
            rμ0=rμ, rμ_max=rμ_cap, grow=grow, tol=tol, verbose=verbose)

        if (!okκ && !okμ)
            break
        end
        rκ, rμ = max(rκ_new, rκ), max(rμ_new, rμ)
    end

    return (rκ > 0 && rμ > 0, κc, μc, rκ, rμ)
end

function cover_rectangle(κ_hi=(0.5 * π)^2, μT=0.0886765;
    nx=50, ny=50, rκ0=1e-9, rμ0=1e-9, verbose=true)

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


# Return the largest validated box sizes (Δκ_max, Δμ_max), plus the radii used
function max_box_sizes(boxes; only_ok::Bool=true)
    rsκ = Float64[]
    rsμ = Float64[]
    for (ok, _κ, _μ, rκ, rμ) in boxes
        if !only_ok || ok
            push!(rsκ, rκ)
            push!(rsμ, rμ)
        end
    end
    isempty(rsκ) && error("No boxes to analyze.")
    rκ_max = maximum(rsκ)
    rμ_max = maximum(rsμ)
    return 2 * rκ_max, 2 * rμ_max, rκ_max, rμ_max
end

function mean_box_sizes(boxes; only_ok::Bool=true)
    rsκ = Float64[]
    rsμ = Float64[]
    for (ok, _κ, _μ, rκ, rμ) in boxes
        if !only_ok || ok
            push!(rsκ, rκ)
            push!(rsμ, rμ)
        end
    end
    isempty(rsκ) && error("No boxes to analyze.")
    rκ_mean = mean(rsκ)
    rμ_mean = mean(rsμ)
    return 2 * rκ_mean, 2 * rμ_mean, rκ_mean, rμ_mean
end


Δκ_max, Δμ_max, rκ_max, rμ_max = max_box_sizes(boxes)
Δκ_mean, Δμ_mean, rκ_mean, rμ_mean = mean_box_sizes(boxes)
println("Mean box width Δκ=", Δκ_mean, " height Δμ=", Δμ_mean)

