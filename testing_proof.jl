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


# function push_kappa_interval(κ_center, μ_center, root;
#     Δκ=1e-3, κmax=(π / 2)^2,
#     r_κ=1e-8, r_μ=1e-8,
#     verbose=false)

#     # first validation
#     ok = validate_f1_root(root, κ_center, μ_center, r_κ, r_μ; verbose)
#     ok || error("initial point not valid")

#     κ_best, root_best = κ_center, root

#     while κ_center < κmax
#         # next center value (float)
#         κ_center_next = κ_center + Δκ

#         # build intervals for rigorous check
#         Iκ = interval(κ_center_next, r_κ; format=:midpoint)
#         Iμ = interval(μ_center, r_μ; format=:midpoint)

#         # compute the root at midpoint values (for seeding)
#         root_next = find_zero(θ -> f₁(θ, mid(Iκ), mid(Iμ)), root)

#         if !isfinite(root_next)
#             verbose && println("nonfinite root at κ=$(mid(Iκ))")
#             break
#         end

#         if validate_f1_root(root_next, mid(Iκ), mid(Iμ), r_κ, r_μ; verbose)
#             κ_center, root = mid(Iκ), root_next
#             κ_best, root_best = κ_center, root
#             verbose && println("✓ valid at κ=$(round(κ_center, sigdigits=6))")
#         else
#             verbose && println("✗ failed at κ=$(round(mid(Iκ), sigdigits=6))")
#             break
#         end
#     end

#     return interval(κ_best, r_κ; format=:midpoint), root_best
# end

# κ0, μ = 0.5, 0.08
# Iκ_final, root_final = push_kappa_interval(κ0, μ, root; Δκ=1e-3, r_κ=1e-8, verbose=true)

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


# """
#     push_forward_mu(initial_root, κc, μc, rκ;
#                     rμ_init=1e-8, μT=0.0886765, rμ_cap=0.1,
#                     grow=1.8, tol=1e-12, max_steps=200, verbose=true)

# Walk forward in μ. At each step:
#   1) maximize Δμ with `maximize_delta_mu` (cap so μ never exceeds μT),
#   2) record the validated box (true, κc, μ_old, rκ, rμ_step),
#   3) move μc ← μ_old + rμ_step,
#   4) find the new root at (κc, μc) using a bracket (falls back to Newton),
#   5) reset rμ ← rμ_init and repeat.

# Returns: `boxes, κc, μc, root`
#   where `boxes :: Vector{Tuple{Bool,Real,Real,Real,Real}}`.
# """
# function push_forward_mu(initial_root, κc, μc, rκ;
#     rμ_init::Real=1e-8,
#     μT::Real=0.0886765,
#     rμ_cap::Real=0.1,        # per-step hard cap (before μT clipping)
#     grow::Real=1.8,
#     tol::Real=1e-12,
#     max_steps::Int=200,
#     verbose::Bool=true)

#     # numeric type hygiene (works for Float64 / BigFloat)
#     T = promote_type(typeof(initial_root), typeof(κc), typeof(μc),
#         typeof(rκ), typeof(rμ_init), typeof(μT), typeof(tol))
#     root = T(initial_root)
#     rμ = T(rμ_init)

#     boxes = Vector{Tuple{Bool,T,T,T,T}}()

#     stepno = 0
#     while μc < μT
#         stepno += 1

#         # don’t allow step to push past μT
#         rμ_max_allowed = min(T(rμ_cap), T(μT - μc))
#         if rμ_max_allowed ≤ T(0)
#             verbose && println("Reached μT = ", μT)
#             break
#         end

#         # maximize μ-radius for this center
#         ok, κc_ret, μc_ret, rκ_ret, rμ_star =
#             maximize_delta_mu(root, κc, μc,
#                 rκ;
#                 rμ0=rμ,
#                 rμ_max=rμ_max_allowed,
#                 grow=grow,
#                 tol=tol,
#                 max_steps=max_steps,
#                 verbose=verbose)

#         if !ok || rμ_star ≤ T(0)
#             verbose && println("Stopping: cannot certify μ-radius at μ = ", μc)
#             break
#         end

#         # save current validated box at the *old* center (plot later)
#         push!(boxes, (true, T(κc), T(μc), T(rκ), T(rμ_star)))

#         # advance μ-center by the certified radius
#         μc_old = μc
#         μc = μc_old + rμ_star
#         verbose && println("→ step ", stepno, ": μ ← ", μc)

#         # find new root at the new center
#         br = brackets_f₁(κc, μc; n=800)
#         if isempty(br)
#             verbose && println("Stopping: could not refine root at μ = ", μc)
#             break
#         end

#         root = find_zero(θ -> f₁(θ, κc, μc), br[1])
#         # reset μ-radius for the next step
#         rμ = T(rμ_init)
#     end

#     return boxes, κc, μc, root
# end



# """
#     sweep_bands_kappa(initial_root, κc0, μ0, rκ;
#         μT = 0.0886765,
#         κ_hi = (0.5*π)^2,
#         rμ_init = 1e-8,
#         rμ_cap  = 0.1,
#         rκ_init = 1e-8,        # starting κ-radius for κ-steps
#         rκ_cap  = 0.1,         # per-step cap when growing κ
#         grow = 1.8,
#         tol  = 1e-12,
#         max_steps = 200,
#         verbose = true)

# Make consecutive μ-bands while increasing κ. For each band:
#   1) Run `push_forward_mu` from (κc, μ_start) with fixed rκ → collect boxes along μ.
#   2) Increase κ by maximizing κ-radius (holding rμ fixed at rμ_init), don’t exceed κ_hi.
#   3) Recompute the root at (κc_new, μ_start). Repeat.

# Returns `(all_boxes, κc, μ_start, last_root)`, where `all_boxes` is a
# Vector{Tuple{Bool,Real,Real,Real,Real}} ready for `plot_validated_boxes`.
# """
# function sweep_bands_kappa(initial_root, κc0, μ0, rκ;
#     μT::Real=0.0886765,
#     κ_hi::Real=(0.5 * π)^2,
#     rμ_init::Real=1e-8,
#     rμ_cap::Real=0.1,
#     rκ_init::Real=1e-8,
#     rκ_cap::Real=0.1,
#     grow::Real=1.8,
#     tol::Real=1e-12,
#     max_steps::Int=200,
#     verbose::Bool=true)

#     T = promote_type(typeof(initial_root), typeof(κc0), typeof(μ0),
#         typeof(rκ), typeof(rμ_init), typeof(μT), typeof(κ_hi), typeof(tol))
#     κc = T(κc0)
#     μ_start = T(μ0)
#     root = T(initial_root)

#     all_boxes = Vector{Tuple{Bool,T,T,T,T}}()

#     band = 0
#     while κc < κ_hi && band < 1000
#         band += 1
#         verbose && println("\n=== Band $band @ κ = $(κc) ===")

#         # 1) μ-push at current κ (collect validated rectangles along μ)
#         boxes, κc_out, μ_end, root_out = push_forward_mu(root, κc, μ_start, rκ;
#             rμ_init=T(rμ_init), μT=T(μT), rμ_cap=T(rμ_cap),
#             grow=T(grow), tol=T(tol), max_steps=max_steps, verbose=verbose)

#         append!(all_boxes, boxes)
#         root = root_out  # last root in this band

#         # stop if we already reached μT inside push_forward_mu (common case)
#         if μ_end < μT && verbose
#             println("Band ended early at μ=$(μ_end).")
#         end

#         # 2) grow κ-radius (hold a small fixed rμ during this κ step)
#         rκ_max_allowed = min(T(rκ_cap), T(κ_hi - κc))
#         if rκ_max_allowed ≤ T(0)
#             verbose && println("Reached κ_hi = ", κ_hi)
#             break
#         end

#         # move κ center
#         κc_old = κc
#         κc = κc_old + rκ
#         verbose && println("→ advance κ: $(κc_old) → $(κc)")

#         # 3) re-find root at the new (κc, μ_start) to seed next band
#         br = brackets_f₁(κc, μ_start; n=200)
#         root = find_zero(θ -> f₁(θ, κc, μc), br[1])

#         if isempty(br)
#             verbose && println("Stopping: could not refine root at (κ,μ) = ($(κc), $(μ_start))")
#             break
#         end

#         okκ, κc_ret, μc_ret, rκ_star, rμ_fixed =
#             maximize_delta_kappa(root, κc, μ_start,
#                 T(rμ_init);      # keep μ-radius tiny during κ move
#                 rκ0=T(rκ_init),
#                 rκ_max=rκ_max_allowed,
#                 grow=T(grow), tol=T(tol), max_steps=max_steps, verbose=verbose)

#                 rκ = rκ_star
#         if !okκ || rκ_star ≤ T(0)
#             verbose && println("Stopping: cannot certify κ-radius at κ = ", κc)
#             break
#         end


#         # next while-iteration: new band at (κc, μ_start)
#     end

#     return all_boxes
# end

# # boxes = sweep_bands_kappa(root, κc, μc, 1e-8;
# #     μT=0.0886765, κ_hi=(0.5 * π)^2, rμ_init=rμ, rμ_cap=0.2,
# #     rκ_init=1e-8, rκ_cap=0.5, grow=1.8, tol=1e-12, verbose=false)

# # # then plot:
# # plot_validated_boxes(boxes; color=:green, alpha=0.5, title="Validated bands in (κ, μ)")


# # plot_validated_boxes(boxes; color=:orange, alpha=0.5, title="Validated regions in (κ, μ)")

# """
#     maximize_box(root, κc, μc; rκ0=1e-8, rμ0=1e-8,
#                  κ_hi=(0.5*π)^2, μT=0.0886765,
#                  grow=1.8, tol=1e-12, itmax=20, verbose=false)

# Alternating maximization of (rκ, rμ) around (κc, μc), stopping when neither radius improves.
# Returns (ok, κc, μc, rκ, rμ).
# """
# function maximize_box(root, κc, μc;
#         rκ0=1e-9, rμ0=1e-9,
#         κ_hi=(0.5*π)^2, μT=0.0886765,
#         grow=1.2, tol=1e-12, itmax=20, verbose=false)

#     rκ, rμ = rκ0, rμ0
#     for t in 1:itmax
#         # cap by remaining distance to boundary
#         rκ_cap = min((κ_hi - κc), κ_hi)  # mild cap; tweak as you wish
#         rμ_cap = min((μT - μc),  μT)

#         # κ step (hold rμ)
#         okκ, _, _, rκ_new, _ = maximize_delta_kappa(root, κc, μc,
#             rμ; rκ0=rκ, rκ_max=rκ_cap, grow=grow, tol=tol, verbose=false)

#         # μ step (hold rκ)
#         okμ, _, _, _, rμ_new = maximize_delta_mu(root, κc, μc,
#             rκ_new; rμ0=rμ, rμ_max=rμ_cap, grow=grow, tol=tol, verbose=false)

#         # stop if nothing improved meaningfully
#         if (!okκ && !okμ) || (abs(rκ_new - rκ) ≤ 1e-3*max(rκ,1e-30) &&
#                               abs(rμ_new - rμ) ≤ 1e-3*max(rμ,1e-30))
#             rκ, rμ = max(rκ,0.0), max(rμ,0.0)
#             ok = rκ > 0 && rμ > 0
#             return (ok, κc, μc, rκ, rμ)
#         end

#         rκ, rμ = rκ_new, rμ_new
#     end
#     (true, κc, μc, rκ, rμ)
# end

# using Roots

# """
#     cover_rectangle(κ_hi=(0.5*π)^2, μT=0.0886765;
#         nx=12, ny=12, rκ0=1e-8, rμ0=1e-8, verbose=true)

# Coarse seeds + local box maximization, returns a vector of (ok, κc, μc, rκ, rμ)
# suitable for plot_validated_boxes().
# """
# function cover_rectangle(κ_hi=(0.5*π)^2, μT=0.0886765;
#         nx=100, ny=100, rκ0=1e-9, rμ0=1e-9, verbose=true)

#     κs = range(0.02, κ_hi; length=nx)
#     μs = range(0.0, μT;  length=ny)

#     boxes = Tuple{Bool,Float64,Float64,Float64,Float64}[]

#     for κc in κs, μc in μs
#         # find a root to seed the local box (skip if bracketing fails)
#         br = brackets_f₁(κc, μc; n=200)
#         if isempty(br)
#             # try a neutral theta guess if you like; otherwise skip
#             verbose && @debug "No bracket at (κ,μ)=($κc,$μc); skipping seed."
#             continue
#         end
#         root = find_zero(θ -> f₁(θ, κc, μc), br[1])

#         ok, κc2, μc2, rκ, rμ = maximize_box(root, κc, μc;
#                                 rκ0=rκ0, rμ0=rμ0, κ_hi=κ_hi, μT=μT, verbose=false)
#         ok && push!(boxes, (ok, κc2, μc2, rκ, rμ))
#     end
#     boxes
# end


# boxes = cover_rectangle()
# plot_validated_boxes(boxes; color=:teal, alpha=0.45, title="Validated cover of (κ, μ)")
