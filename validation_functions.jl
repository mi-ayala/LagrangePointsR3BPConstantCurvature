ispos(I) = inf(I) > 0
isneg(I) = sup(I) < 0
nonzero(I) = ispos(I) || isneg(I)

# using IntervalArithmetic

# # helpers (only if not already in your code)
# ispos(I) = isfinite(I) && inf(I) > 0
# isneg(I) = isfinite(I) && sup(I) < 0
# nonzero(I) = ispos(I) || isneg(I)

"""
    certify_deriv_sign_centered(θc, κ, μ, rκ, rμ; min_width=1e-12, verbose=false)

Centered at θc, shrink r by halves within Iθ = [0, θ1_max] until we can certify
f₁θ > 0 or < 0 on Ic = [θc-r, θc+r]. First try a direct call f₁θ(Ic, Iκ, Iμ);
only if that fails, call eval_f₁θ_subdiv(Ic, Iκ, Iμ).
"""
function certify_deriv_sign_centered(θc, κ, μ, rκ, rμ; min_width=1e-12, verbose=false)
    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    Iθ = interval(0, θ1_max(Iκ, Iμ))
    inf(Iθ) < θc < sup(Iθ) || return (:unknown, interval(θc), Interval(0))

    rmax = min(θc - inf(Iθ), sup(Iθ) - θc)
    r = rmax

    for _ in 0:60
        Ic = interval(θc, r; format=:midpoint)

        # 1) try direct evaluation first
        dI = f₁θ(Ic, Iκ, Iμ)
        if ispos(dI)
            verbose && println("f₁θ > 0 on ", Ic, " (direct)")
            return (:pos, Ic, dI)
        elseif isneg(dI)
            verbose && println("f₁θ < 0 on ", Ic, " (direct)")
            return (:neg, Ic, dI)
        end

        # 2) fallback: subdivided evaluation
        dIs = eval_f₁θ_subdiv(Ic, Iκ, Iμ)
        if ispos(dIs)
            verbose && println("f₁θ > 0 on ", Ic, " (subdiv)")
            return (:pos, Ic, dIs)
        elseif isneg(dIs)
            verbose && println("f₁θ < 0 on ", Ic, " (subdiv)")
            return (:neg, Ic, dIs)
        end

        if 2r <= min_width
            verbose && println("inconclusive at width≈", 2r, "  dI=", dIs)
            return (:unknown, Ic, dIs)
        end
        r *= 0.5
    end

    Ic = interval(θc, r; format=:midpoint)
    return (:unknown, Ic, f₁θ(Ic, Iκ, Iμ))
end


"""
    validate_complements(I, J, κ, μ, rκ, rμ; verbose=false)

If J ≠ I, certify f₁ ≠ 0 on left [a,c] and right [d,b].
First try direct f₁ on each complement; only if inconclusive, fallback to eval_f₁_subdiv.
"""
function validate_complements(I, J, κ, μ, rκ, rμ; verbose=false)
    a, b = inf(I), sup(I)
    c, d = inf(J), sup(J)
    a ≤ c ≤ d ≤ b || error("J must be contained in I")

    if a == c && b == d
        verbose && println("J equals I — nothing to check.")
        return true
    end

    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    # Left complement
    ok_left = true
    if a < c
        left = interval(a, c)
        # direct try
        valL = f₁(left, Iκ, Iμ)
        if !nonzero(valL)
            # fallback: subdivided
            valL = eval_f₁_subdiv(left, Iκ, Iμ)
        end
        ok_left = nonzero(valL)
        verbose && println("Left [a,c] = $left → f₁(left) = $valL  → nonzero? ", ok_left)
    end

    # Right complement
    ok_right = true
    if d < b
        right = interval(d, b)
        # direct try
        valR = f₁(right, Iκ, Iμ)
        if !nonzero(valR)
            # fallback: subdivided
            valR = eval_f₁_subdiv(right, Iκ, Iμ)
        end
        ok_right = nonzero(valR)
        verbose && println("Right [d,b] = $right → f₁(right) = $valR  → nonzero? ", ok_right)
    end

    result = ok_left && ok_right
    verbose && println("Complements nonzero? ", result)
    return result
end


using IntervalArithmetic

# ---- f₁θ version ----
function eval_f₁θ_subdiv(Iθ::Interval, Iκ::Interval, Iμ::Interval; ntiles::Int=16)
    a, b = inf(Iθ), sup(Iθ)
    step = (b - a) / ntiles
    H = f₁θ(interval(a, a + step), Iκ, Iμ)
    for k in 1:ntiles-1
        J = interval(a + k*step, a + (k+1)*step)
        H = hull(H, f₁θ(J, Iκ, Iμ))
    end
    return H
end

# ---- f₁ version ----
function eval_f₁_subdiv(Iθ::Interval, Iκ::Interval, Iμ::Interval; ntiles::Int=16)
    a, b = inf(Iθ), sup(Iθ)
    step = (b - a) / ntiles
    H = f₁(interval(a, a + step), Iκ, Iμ)
    for k in 1:ntiles-1
        J = interval(a + k*step, a + (k+1)*step)
        H = hull(H, f₁(J, Iκ, Iμ))
    end
    return H
end






# function certify_deriv_sign_centered(θc, κ, μ, rκ, rμ; min_width=1e-12, verbose=false)

#     Iκ = interval(κ, rκ; format=:midpoint)
#     Iμ = interval(μ, rμ; format=:midpoint)

#     Iθ = interval(0, θ1_max(Iκ, Iμ))
#     inf(Iθ) < θc < sup(Iθ) || return (:unknown, interval(θc), Interval(0))  # θc must be inside



#     # maximum admissible radius that keeps [θc - r, θc + r] ⊆ Iθ
#     rmax = min(θc - inf(Iθ), sup(Iθ) - θc)
#     r = rmax

#     # shrink-by-halves until strict sign is certified or too small
#     for _ in 0:60
#         Ic = interval(θc, r; format=:midpoint)
#         dI = eval_f₁θ_subdiv(Ic, Iκ, Iμ)
#         if ispos(dI)
#             verbose && println("f₁θ > 0 on ", Ic)
#             return (:pos, Ic, dI)
#         elseif isneg(dI)
#             verbose && println("f₁θ < 0 on ", Ic)
#             return (:neg, Ic, dI)
#         end
#         if 2r <= min_width
#             verbose && println("inconclusive at width≈", 2r, "  dI=", dI)
#             return (:unknown, Ic, dI)
#         end
#         r *= 0.5
#     end
#     Ic = interval(θc, r; format=:midpoint)
#     return (:unknown, Ic, f₁θ(Ic, Iκ, Iμ))
# end


# """
#     validate_complements(I, J, κ, μ, rκ, rμ; verbose=false)

# Checks if J == I. If not, verifies f₁(θ, κ, μ) ≠ 0 on both complements:
# left [a,c] and right [d,b].

# Returns `true` if:
#   • J == I, or
#   • f₁ is certified nonzero on both complement intervals.

# Otherwise returns `false`.
# """
# function validate_complements(I, J, κ, μ, rκ, rμ; verbose=false)
#     a, b = inf(I), sup(I)
#     c, d = inf(J), sup(J)
#     a ≤ c ≤ d ≤ b || error("J must be contained in I")

#     # If exactly equal → done
#     if a == c && b == d
#         verbose && println("J equals I — nothing to check.")
#         return true
#     end

#     Iκ = interval(κ, rκ; format=:midpoint)
#     Iμ = interval(μ, rμ; format=:midpoint)

#     # Left and right complements (skip if empty)
#     left = (a < c) ? interval(a, c) : nothing
#     right = (d < b) ? interval(d, b) : nothing

#     ok_left = true
#     ok_right = true

#     if left !== nothing
#         valL = eval_f₁_subdiv(left, Iκ, Iμ)
#         ok_left = nonzero(valL)
#         verbose && println("Left [a,c] = $left → f₁(left) = $valL")
#     end

#     if right !== nothing
#         valR = eval_f₁_subdiv(right, Iκ, Iμ)
#         ok_right = nonzero(valR)
#         verbose && println("Right [d,b] = $right → f₁(right) = $valR")
#     end

#     result = ok_left && ok_right
#     verbose && println("Complements nonzero? ", result)
#     return result
# end


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

    leftθ = (inf(Iθ) + root) / 2
    rightθ = (root + sup(Iθ)) / 2
    fL = f₁(interval(leftθ), Iκ, Iμ)
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


