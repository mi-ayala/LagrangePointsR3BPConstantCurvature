
ispos(I) = inf(I) > 0
isneg(I) = sup(I) < 0
nonzero(I) = ispos(I) || isneg(I)


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

function eval_f₁θ_subdiv(Iθ::Interval, Iκ::Interval, Iμ::Interval; ntiles::Int=16)
    a, b = inf(Iθ), sup(Iθ)
    step = (b - a) / ntiles
    H = f₁θ(interval(a, a + step), Iκ, Iμ)
    for k in 1:ntiles-1
        J = interval(a + k * step, a + (k + 1) * step)
        H = hull(H, f₁θ(J, Iκ, Iμ))
    end
    return H
end

function eval_f₁_subdiv(Iθ::Interval, Iκ::Interval, Iμ::Interval; ntiles::Int=16)
    a, b = inf(Iθ), sup(Iθ)
    step = (b - a) / ntiles
    H = f₁(interval(a, a + step), Iκ, Iμ)
    for k in 1:ntiles-1
        J = interval(a + k * step, a + (k + 1) * step)
        H = hull(H, f₁(J, Iκ, Iμ))
    end
    return H
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

################################
#### Validation function for f4
################################

"""
    certify_deriv_sign_centered(θc, κ, μ, rκ, rμ; min_width=1e-12, verbose=false)

Centered at θc, shrink r by halves within Iθ = [0, θ1_max] until we can certify
f₁θ > 0 or < 0 on Ic = [θc-r, θc+r]. First try a direct call f₁θ(Ic, Iκ, Iμ);
only if that fails, call eval_f₁θ_subdiv(Ic, Iκ, Iμ).
"""
function certify_deriv_sign_centered_f4(θc, κ, μ, rκ, rμ; min_width=1e-12, verbose=false)
    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    Iθ = interval(inf(θ4_min(Iκ, Iμ)), sup(θ4_max(Iκ, Iμ)))
    inf(Iθ) < θc < sup(Iθ) || return (:unknown, interval(θc), Interval(0))

    rmax = min(θc - inf(Iθ), sup(Iθ) - θc)
    r = rmax

    for _ in 0:60
        Ic = interval(θc, r; format=:midpoint)

        # 1) try direct evaluation first
        dI = f₄θ(Ic, Iκ, Iμ)
        if ispos(dI)
            verbose && println("f₄θ > 0 on ", Ic, " (direct)")
            return (:pos, Ic, dI)
        elseif isneg(dI)
            verbose && println("f₄θ < 0 on ", Ic, " (direct)")
            return (:neg, Ic, dI)
        end

        # 2) fallback: subdivided evaluation
        dIs = eval_f₄θ_subdiv(Ic, Iκ, Iμ)
        if ispos(dIs)
            verbose && println("f₄θ > 0 on ", Ic, " (subdiv)")
            return (:pos, Ic, dIs)
        elseif isneg(dIs)
            verbose && println("f₄θ < 0 on ", Ic, " (subdiv)")
            return (:neg, Ic, dIs)
        end

        if 2r <= min_width
            verbose && println("inconclusive at width≈", 2r, "  dI=", dIs)
            return (:unknown, Ic, dIs)
        end
        r *= 0.5
    end

    Ic = interval(θc, r; format=:midpoint)
    return (:unknown, Ic, f₄θ(Ic, Iκ, Iμ))
end


"""
    validate_complements(I, J, κ, μ, rκ, rμ; verbose=false)

If J ≠ I, certify f₁ ≠ 0 on left [a,c] and right [d,b].
First try direct f₁ on each complement; only if inconclusive, fallback to eval_f₁_subdiv.
"""
function validate_complements_f4(I, J, κ, μ, rκ, rμ; verbose=false)
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
        valL = f₄(left, Iκ, Iμ)
        if !nonzero(valL)
            # fallback: subdivided
            valL = eval_f₄_subdiv(left, Iκ, Iμ)
        end
        ok_left = nonzero(valL)
        verbose && println("Left [a,c] = $left → f₄(left) = $valL  → nonzero? ", ok_left)
    end

    # Right complement
    ok_right = true
    if d < b
        right = interval(d, b)
        # direct try
        valR = f₄(right, Iκ, Iμ)
        if !nonzero(valR)
            # fallback: subdivided
            valR = eval_f₄_subdiv(right, Iκ, Iμ)
        end
        ok_right = nonzero(valR)
        verbose && println("Right [d,b] = $right → f₄(right) = $valR  → nonzero? ", ok_right)
    end

    result = ok_left && ok_right
    verbose && println("Complements nonzero? ", result)
    return result
end

function eval_f₄θ_subdiv(Iθ::Interval, Iκ::Interval, Iμ::Interval; ntiles::Int=16)
    a, b = inf(Iθ), sup(Iθ)
    step = (b - a) / ntiles
    H = f₄θ(interval(a, a + step), Iκ, Iμ)
    for k in 1:ntiles-1
        J = interval(a + k * step, a + (k + 1) * step)
        H = hull(H, f₄θ(J, Iκ, Iμ))
    end
    return H
end

function eval_f₄_subdiv(Iθ::Interval, Iκ::Interval, Iμ::Interval; ntiles::Int=16)
    a, b = inf(Iθ), sup(Iθ)
    step = (b - a) / ntiles
    H = f₄(interval(a, a + step), Iκ, Iμ)
    for k in 1:ntiles-1
        J = interval(a + k * step, a + (k + 1) * step)
        H = hull(H, f₄(J, Iκ, Iμ))
    end
    return H
end

function validate_f4_root(root, κ, μ, r_κ, r_μ; verbose=false)
    # parameter intervals
    Iκ = interval(κ, r_κ; format=:midpoint)
    Iμ = interval(μ, r_μ; format=:midpoint)

    # theta interval [θ4_min, θ4_max] (conservative upper bound)
    a = θ4_min(Iκ, Iμ)
    b = θ4_max(Iκ, Iμ)
    Iθ = interval(inf(a), sup(b))

    # guard: root ∈ Iθ
    if !(inf(Iθ) ≤ root ≤ sup(Iθ))
        verbose && println("root=$root not inside θ-interval $Iθ")
        return false
    end

    leftθ = (inf(Iθ) + root) / 2
    rightθ = (root + sup(Iθ)) / 2
    fL = eval_f₄_subdiv(interval(leftθ), Iκ, Iμ)
    fR = eval_f₄_subdiv(interval(rightθ), Iκ, Iμ)
    left_neg_right_pos = isneg(fL) && ispos(fR)
    left_pos_right_neg = ispos(fL) && isneg(fR)
    exists = left_neg_right_pos || left_pos_right_neg
    if !exists
        verbose && println("Existence not certified by sign change: fL=$fL, fR=$fR")
        return false
    end

    # --- replace old derivative check with your two functions ---
    # 1) certify strict sign of f₄θ on a centered sub-interval J around root
    status, J, dI = certify_deriv_sign_centered_f4(root, κ, μ, r_κ, r_μ)
    if status === :unknown
        verbose && println("Could not certify strict sign of f₄θ near root. dI = $dI")
        return false
    end

    # 2) prove f₄ ≠ 0 on the complements Iθ \ J
    uniq = validate_complements_f4(Iθ, J, κ, μ, r_κ, r_μ; verbose=verbose)

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


##################################
### Stability validation functions
##################################

# [0,0] with the same numeric type as I
@inline zero_interval_like(I::Interval) = Interval(zero(inf(I)))

# Robust sign tagging with no `==` on intervals and no `in` fallback
@inline function sign_tag(I::Interval)
    if ispos(I)
        return :pos
    elseif isneg(I)
        return :neg
    elseif isequal_interval(I, zero_interval_like(I))
        return :zero
    elseif (inf(I) ≤ 0) && (0 ≤ sup(I))
        return :mixed   # contains 0 but not thin zero
    else
        return :unknown
    end
end

"""
    classify_equilibrium(Iθ, Iκ, Iμ; verbose=false)

Implements the 3-step test on intervals:

1) If det(θ) < 0  → :saddle_hyperbolic
2) else if λ₂(θ) > 0 → :lyapunov_stable
3) else if Δ(θ) > 0 → :elliptic
   else if Δ(θ) < 0 → :linearly_unstable
   else → :unknown

Returns `(ok::Bool, label::Symbol, Idet, Ilam2, IΔ)`.
`ok == true` iff the label is certified (not :unknown).
"""
function classify_equilibrium(Iθ, Iκ, Iμ; verbose::Bool=false)
    Idet = detλ(Iθ, Iκ, Iμ)
    tgdet = sign_tag(Idet)

    if verbose
        println("det(θ) ∈ ", Idet, "  ⇒ ", tgdet)
    end

    if tgdet === :neg
        return true, :saddle_hyperbolic, Idet, IntervalArithmetic.interval(0), IntervalArithmetic.interval(0)
    end

    Ilam2 = λ₂(Iθ, Iκ, Iμ)
    tgl2 = sign_tag(Ilam2)
    if verbose
        println("λ₂(θ) ∈ ", Ilam2, "  ⇒ ", tgl2)
    end

    if tgl2 === :pos
        return true, :lyapunov_stable, Idet, Ilam2, IntervalArithmetic.interval(0)
    end

    IΔ = Δ(Iθ, Iκ, Iμ)
    tgΔ = sign_tag(IΔ)
    if verbose
        println("Δ(θ) ∈ ", IΔ, "  ⇒ ", tgΔ)
    end

    if tgΔ === :pos
        return true, :elliptic, Idet, Ilam2, IΔ
    else
        return true, :linearly_unstable, Idet, Ilam2, IΔ
    end
end


"""
    refine_root_interval(root, κ, μ, rκ, rμ, J;
                         steps=10, verbose=false)

Given:
  - root  :: Real    (numerical root)
  - κ, μ  :: Reals   (box center)
  - rκ,rμ :: Reals   (parameter radii)
  - J     :: Interval  (certified existence+uniqueness interval for θ)

Return a smaller interval Jref ⊂ J that still rigorously contains the
unique root, obtained by steps sign-change checks around `root`.
"""
function refine_root_interval(root, κ, μ, rκ, rμ, J::Interval;
    steps::Int=10, verbose::Bool=false)

    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    left = inf(J)
    right = sup(J)

    # ensure root is inside J (it should be, if earlier code is correct)
    if !(left ≤ root ≤ right)
        verbose && println("root=$root not inside certified J=$J; no refinement.")
        return J
    end

    for k in 1:steps
        width = right - left
        width ≤ 0 && break

        # propose a symmetric shrink around root
        δ = width / 4          # shrink factor; can tune
        new_left = root - δ
        new_right = root + δ

        # keep new interval inside original J
        if new_left <= left || new_right >= right
            verbose && println("Step $k: proposed interval not inside J; stop.")
            break
        end

        fL = f₁(interval(new_left), Iκ, Iμ)
        fR = f₁(interval(new_right), Iκ, Iμ)

        left_neg_right_pos = isneg(fL) && ispos(fR)
        left_pos_right_neg = ispos(fL) && isneg(fR)
        has_sign_change = left_neg_right_pos || left_pos_right_neg

        if has_sign_change
            verbose && println("Step $k: shrink to [$new_left, $new_right], fL=$fL, fR=$fR")
            left, right = new_left, new_right
        else
            verbose && println("Step $k: cannot shrink further (no clean sign change).")
            break
        end
    end

    Jref = interval(left, right)
    verbose && println("Refined J: ", Jref)
    return Jref
end


"""
    validate_f1_root_with_stability(root, κ, μ, rκ, rμ;
                                    refine_steps::Int=8,
                                    verbose::Bool=false)

1. Prove existence of a root of f₁(θ, κ, μ) in [0, θ₁_max] via sign change.
2. Prove uniqueness via derivative sign on a centered subinterval and
   non-vanishing of f₁ on the complements.
3. Refine the certified θ-interval around `root` using `refine_root_interval`.
4. Classify stability on the refined θ-interval using `classify_equilibrium`.

Returns:
    (ok::Bool, label::Symbol, J_ref, Idet, Ilam2, IΔ)

- `ok` is true iff existence+uniqueness AND a stability classification
  decision (possibly :unknown) were successfully performed.
- `label ∈ {:saddle_hyperbolic, :lyapunov_stable, :elliptic,
            :linearly_unstable, :unknown}`
- `J_ref` is the refined existence interval for θ.
- `Idet, Ilam2, IΔ` are the interval enclosures used in the classification.
"""
function validate_f1_root_with_stability(root, κ, μ, rκ, rμ;
    refine_steps::Int=8,
    verbose::Bool=false)

    # --- 1) parameter intervals ---
    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    # θ-interval [0, θ₁_max] (conservative)
    a = zero(κ)
    b = θ1_max(Iκ, Iμ)
    Iθ = interval(a, sup(b))

    # guard: root ∈ Iθ
    if !(inf(Iθ) ≤ root ≤ sup(Iθ))
        verbose && println("root=$root not inside θ-interval $Iθ")
        return false, :unknown, Iθ, Interval(0), Interval(0), Interval(0)
    end

    # --- 2) existence by sign change at midpoints ---
    leftθ = (inf(Iθ) + root) / 2
    rightθ = (root + sup(Iθ)) / 2
    fL = f₁(interval(leftθ), Iκ, Iμ)
    fR = f₁(interval(rightθ), Iκ, Iμ)

    left_neg_right_pos = isneg(fL) && ispos(fR)
    left_pos_right_neg = ispos(fL) && isneg(fR)
    exists = left_neg_right_pos || left_pos_right_neg

    if !exists
        verbose && println("Existence not certified by sign change: fL=$fL, fR=$fR")
        return false, :unknown, Iθ, Interval(0), Interval(0), Interval(0)
    end

    # --- 3) uniqueness: derivative sign + complements ---
    status, J, dI = certify_deriv_sign_centered(root, κ, μ, rκ, rμ; verbose=false)
    if status === :unknown
        verbose && println("Could not certify strict sign of f₁θ near root. dI = $dI")
        return false, :unknown, Iθ, Interval(0), Interval(0), Interval(0)
    end

    uniq = validate_complements(Iθ, J, κ, μ, rκ, rμ; verbose=verbose)
    if !uniq
        verbose && println("Uniqueness via complements failed.")
        return false, :unknown, J, Interval(0), Interval(0), Interval(0)
    end

    verbose && begin
        println("θ-interval Iθ     = $Iθ")
        println("leftθ=$leftθ  fL  = $fL  (neg? ", isneg(fL), ", pos? ", ispos(fL), ")")
        println("rightθ=$rightθ fR = $fR  (neg? ", isneg(fR), ", pos? ", ispos(fR), ")")
        println("certified J       = $J (status = $status)")
        println("uniqueness via complements? ", uniq)
    end

    # --- 4) refine the existence interval around root (cheap) ---
    J_ref = refine_root_interval(root, κ, μ, rκ, rμ, J;
        steps=refine_steps, verbose=verbose)

    # --- 5) classify stability on J_ref ---
    okc, label, Idet, Ilam2, IΔ = classify_equilibrium(J_ref, Iκ, Iμ; verbose=verbose)

    # ok = existence+uniqueness succeeded, classification ran (even if label=:unknown)
    ok_total = uniq && okc

    return ok_total, label, J_ref, Idet, Ilam2, IΔ
end

function refine_root_interval_f4(root, κ, μ, rκ, rμ, J::Interval;
    steps::Int=10, verbose::Bool=false)

    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    left = inf(J)
    right = sup(J)

    # ensure root is inside J (it should be, if earlier code is correct)
    if !(left ≤ root ≤ right)
        verbose && println("root=$root not inside certified J=$J; no refinement.")
        return J
    end

    for k in 1:steps
        width = right - left
        width ≤ 0 && break

        # propose a symmetric shrink around root
        δ = width / 4          # shrink factor; can tune
        new_left = root - δ
        new_right = root + δ

        # keep new interval inside original J
        if new_left <= left || new_right >= right
            verbose && println("Step $k: proposed interval not inside J; stop.")
            break
        end

        fL = f₄(interval(new_left), Iκ, Iμ)
        fR = f₄(interval(new_right), Iκ, Iμ)

        left_neg_right_pos = isneg(fL) && ispos(fR)
        left_pos_right_neg = ispos(fL) && isneg(fR)
        has_sign_change = left_neg_right_pos || left_pos_right_neg

        if has_sign_change
            verbose && println("Step $k: shrink to [$new_left, $new_right], fL=$fL, fR=$fR")
            left, right = new_left, new_right
        else
            verbose && println("Step $k: cannot shrink further (no clean sign change).")
            break
        end
    end

    Jref = interval(left, right)
    verbose && println("Refined J: ", Jref)
    return Jref
end