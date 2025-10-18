ispos(I) = inf(I) > 0
isneg(I) =  sup(I) < 0
nonzero(I) = ispos(I) || isneg(I)

# ispos(I) = inf(I) > 0
# isneg(I) = sup(I) < 0
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


function check_uniqueness(κ, μ, r_κ, r_μ; verbose=false)
    Iκ = interval(κ, r_κ; format=:midpoint)
    Iμ = interval(μ, r_μ; format=:midpoint)
    Iθ = interval(0, θ1_max(Iκ, Iμ))

    status, Iθ_cert, dI_cert = certify_deriv_sign_centered(θc, κ, μ, rκ, rμ; verbose=false)

    ok = validate_complements(Iθ, Iθ_cert, κ, μ, rκ, rμ; verbose=true)

    return ok
end

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