ispos(I) = inf(I) > 0
isneg(I) =  sup(I) < 0
nonzero(I) = ispos(I) || isneg(I)


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

    # theta interval [0, θ1_max] (conservative upper bound)
    a = zero(κ)
    b = θ1_max(Iκ, Iμ)
    Iθ = interval(a, sup(b))

    # guard: root ∈ Iθ
    if !(inf(Iθ) ≤ root ≤ sup(Iθ))
        verbose && println("root=$root not inside θ-interval $Iθ")
        return false
    end

    leftθ  = (inf(Iθ) + root) / 2
    rightθ = (root + sup(Iθ)) / 2
    fL = f₁(interval(leftθ),  Iκ, Iμ)
    fR = f₁(interval(rightθ), Iκ, Iμ)
    left_neg_right_pos = isneg(fL) && ispos(fR)
    left_pos_right_neg = ispos(fL) && isneg(fR)
    exists = left_neg_right_pos || left_pos_right_neg
    if !exists
        verbose && println("Existence not certified by sign change: fL=$fL, fR=$fR")
        return false
    end

    # --- replace old derivative check with your two functions ---
    # 1) certify strict sign of f₁θ on a centered sub-interval J around root
    status, J, dI = certify_deriv_sign_centered(root, κ, μ, r_κ, r_μ)
    if status === :unknown
        verbose && println("Could not certify strict sign of f₁θ near root. dI = $dI")
        return false
    end

    # 2) prove f₁ ≠ 0 on the complements Iθ \ J
    uniq = validate_complements(Iθ, J, κ, μ, r_κ, r_μ; verbose=verbose)

    if verbose
        println("θ-interval Iθ     = $Iθ")
        println("leftθ=$leftθ  fL  = $fL  (neg? ", isneg(fL), ", pos? ", ispos(fL), ")")
        println("rightθ=$rightθ fR = $fR  (neg? ", isneg(fR), ", pos? ", ispos(fR), ")")
        println("certified J       = $J (status = $status)")
        println("uniqueness via complements? ", uniq)
    end

    # return (uniq, κ, μ, r_κ, r_μ)
    return uniq
end


