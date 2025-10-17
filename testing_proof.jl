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

# f₁(θ, κ, μ)
# f₁θ(θ, κ, μ)
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


function θ1_max(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    y = tan(r * q₂(κ, μ))     # numerator
    x = one(μ) + μ            # denominator (kept as same type)
    return atan(y, x)
end



f₁θ(interval(0, θ1_max(interval(κ), interval(μ))), interval(κ), interval(μ))


# # bounds for the box
# μT = 0.0886765
# κ_lo, κ_hi = 0.0, (0.5 * pi)^2
# μ_lo, μ_hi = 0.0, μT



# ### A function that validates a root


# # helpers for interval sign logic
# ispos(I) = isfinite(I) && inf(I) > 0
# isneg(I) = isfinite(I) && sup(I) < 0


# helpers for interval sign logic
ispos(I) = inf(I) > 0
isneg(I) = sup(I) < 0
"""
    validate_f1_root(root, κ, μ, r_κ, r_μ; verbose=false, shrink=1e-12)

Validate that `root` is a simple root of f₁(·; κ, μ) in θ ∈ [0, θ1_max(κ, μ)].
`κ` and `μ` are centered to intervals via radii `r_κ`, `r_μ` (midpoint format).

Checks:
  - f₁(left) and f₁(right) have opposite signs (using thin intervals at midpoints)
  - the derivative f₁θ is monotone over the full θ-interval (strictly >0 or <0)

Returns `true`/`false`. If `verbose=true`, prints diagnostic info.
"""
function validate_f1_root(root, κ, μ, r_κ, r_μ; verbose=false)
    # parameter intervals
    Iκ = interval(κ, r_κ; format=:midpoint)
    Iμ = interval(μ, r_μ; format=:midpoint)

    # theta interval [0, b], slightly shrunk to avoid boundary singularities
    a = zero(κ)
    b = θ1_max(Iκ, Iμ)
    A = interval(a, sup(b))  # conservative

    # guard: root inside [A]?
    if !(inf(A) ≤ root ≤ sup(A))
        verbose && println("root=$root not inside θ-interval $A")
        return false
    end

    # midpoints on each side of root
    leftθ = (inf(A) + root) / 2
    rightθ = (root + sup(A)) / 2

    # thin-interval evaluations
    fL = f₁(interval(leftθ), Iκ, Iμ)
    fR = f₁(interval(rightθ), Iκ, Iμ)

    # derivative over the full interval
    dI = f₁θ(A, Iκ, Iμ)

    left_neg_right_pos = isneg(fL) && ispos(fR)
    left_pos_right_neg = ispos(fL) && isneg(fR)
    deriv_pos = ispos(dI)
    deriv_neg = isneg(dI)

    ok = (left_neg_right_pos && deriv_pos) || (left_pos_right_neg && deriv_neg)

    if verbose
        println("θ-interval A       = $A")
        println("leftθ=$leftθ  fL   = $fL  (neg? ", isneg(fL), ", pos? ", ispos(fL), ")")
        println("rightθ=$rightθ fR  = $fR  (neg? ", isneg(fR), ", pos? ", ispos(fR), ")")
        println("f₁θ over A         = $dI  (deriv>0? ", deriv_pos, ", deriv<0? ", deriv_neg, ")")
        println("VALID? ", ok)
    end

    return ok
end

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


ispos(I) = isfinite(I) && inf(I) > 0
isneg(I) = isfinite(I) && sup(I) < 0


# Certify sign of f₁θ on a centered sub-interval around θc, clipped to stay inside Iθ
function certify_deriv_sign_centered(θc, κ, μ, rκ, rμ; min_width=1e-12, verbose=false)

    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    Iθ = interval(0, θ1_max(Iκ, Iμ))
    inf(Iθ) < θc < sup(Iθ) || return (:unknown, interval(θc), Interval(0))  # θc must be inside



    # maximum admissible radius that keeps [θc - r, θc + r] ⊆ Iθ
    rmax = min(θc - inf(Iθ), sup(Iθ) - θc)
    r = rmax

    # shrink-by-halves until strict sign is certified or too small
    for _ in 0:60
        Ic = interval(θc, r; format=:midpoint)
        dI = f₁θ(Ic, Iκ, Iμ)
        if ispos(dI)
            verbose && println("f₁θ > 0 on ", Ic)
            return (:pos, Ic, dI)
        elseif isneg(dI)
            verbose && println("f₁θ < 0 on ", Ic)
            return (:neg, Ic, dI)
        end
        if 2r <= min_width
            verbose && println("inconclusive at width≈", 2r, "  dI=", dI)
            return (:unknown, Ic, dI)
        end
        r *= 0.5
    end
    Ic = interval(θc, r; format=:midpoint)
    return (:unknown, Ic, f₁θ(Ic, Iκ, Iμ))
end



θc = root
rθ = interval(0, θ1_max(interval(κ), interval(μ)))
κ = 0.5
μ = 0.02
rκ = 1e-8
rμ = 1e-8

status, Iθ_cert, dI_cert = certify_deriv_sign_centered(θc, κ, μ, rκ, rμ; verbose=true)



ispos(I) = isfinite(I) && inf(I) > 0
isneg(I) = isfinite(I) && sup(I) < 0
nonzero(I) = ispos(I) || isneg(I)

"""
    validate_complements(I, J, κ, μ, rκ, rμ; verbose=false)

Checks if J == I. If not, verifies f₁(θ, κ, μ) ≠ 0 on both complements:
left [a,c] and right [d,b].

Returns `true` if:
  • J == I, or
  • f₁ is certified nonzero on both complement intervals.

Otherwise returns `false`.
"""
function validate_complements(I, J, κ, μ, rκ, rμ; verbose=false)
    a, b = inf(I), sup(I)
    c, d = inf(J), sup(J)
    a ≤ c ≤ d ≤ b || error("J must be contained in I")

    # If exactly equal → done
    if a == c && b == d
        verbose && println("J equals I — nothing to check.")
        return true
    end

    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    # Left and right complements (skip if empty)
    left = (a < c) ? interval(a, c) : nothing
    right = (d < b) ? interval(d, b) : nothing

    ok_left = true
    ok_right = true

    if left !== nothing
        valL = f₁(left, Iκ, Iμ)
        ok_left = nonzero(valL)
        verbose && println("Left [a,c] = $left → f₁(left) = $valL")
    end

    if right !== nothing
        valR = f₁(right, Iκ, Iμ)
        ok_right = nonzero(valR)
        verbose && println("Right [d,b] = $right → f₁(right) = $valR")
    end

    result = ok_left && ok_right
    verbose && println("Complements nonzero? ", result)
    return result
end


Iκ = interval(κ, rκ; format=:midpoint)
Iμ = interval(μ, rμ; format=:midpoint)

Iθ = interval(0, θ1_max(Iκ, Iμ))

Jθ = Iθ_cert
κ, μ = 0.5, 0.02
rκ, rμ = 1e-8, 1e-8

ok = validate_complements(Iθ, Jθ, κ, μ, rκ, rμ; verbose=true)
