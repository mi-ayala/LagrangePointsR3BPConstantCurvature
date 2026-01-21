using IntervalArithmetic


ispos(I) = inf(I) > 0
isneg(I) = sup(I) < 0
nonzero(I) = ispos(I) || isneg(I)


"""
    eval_subdiv(f, Iθ::Interval, Iκ::Interval, Iμ::Interval; ntiles::Int=16)

Evaluate an interval function `f(Iθ, Iκ, Iμ)` by subdividing `Iθ` into
`ntiles` equal subintervals and taking the hull of all evaluations.

This is the simplest, non-adaptive scheme:
    [a,b] → [a, a+h], [a+h, a+2h], …, [a+(n-1)h, b]

Returns the hull of f on all tiles.
"""
function eval_subdiv(f::Function,
    Iθ::Interval,
    Iκ::Interval,
    Iμ::Interval; ntiles::Int=16)

    a, b = inf(Iθ), sup(Iθ)
    if ntiles <= 0 || a == b
        return f(Iθ, Iκ, Iμ)
    end

    h = (b - a) / ntiles

    # first tile
    left = a
    right = a + h
    H = f(interval(left, right), Iκ, Iμ)

    # remaining tiles
    for k in 1:ntiles-1
        left = a + k * h
        right = (k == ntiles - 1) ? b : a + (k + 1) * h  # last tile closes at b
        J = interval(left, right)
        H = hull(H, f(J, Iκ, Iμ))
    end

    return H
end

"""
    certify_sign_interval_centered(
        f, θc, Iθ, Iκ, Iμ;
        ntiles=16,
        min_width=1e-12,
        max_halvings=60,
        verbose=false,
    )

Given:
  • f(Iθloc, Iκ, Iμ)         – interval evaluation of a scalar quantity
  • θc                       – center point (typically an approximate root)
  • Iθ                       – allowed θ-interval
  • Iκ, Iμ                   – parameter intervals

Try symmetric intervals Ic = [θc-r, θc+r] ⊆ Iθ, starting from the largest
possible radius r and halving r until we can certify that `f` has strict
sign on Ic:

  • `ispos(FI)` ⇒ f > 0 on Ic
  • `isneg(FI)` ⇒ f < 0 on Ic

For each Ic:
  1) First try a direct call `f(Ic, Iκ, Iμ)`.
  2) If that interval contains 0, refine via `eval_subdiv(f, Ic, Iκ, Iμ; ntiles)`.

Stops when:
  • strict sign is certified, or
  • width(Ic) ≲ min_width (inconclusive).

Returns:
    (status::Symbol, Ic::Interval, FI::Interval)

where:

  • status = :pos  → f > 0 on Ic
  • status = :neg  → f < 0 on Ic
  • status = :unknown → no strict sign certified down to width ≈ min_width
"""
function certify_sign_interval_centered(
    f::Function,
    θc,
    Iθ::Interval,
    Iκ::Interval,
    Iμ::Interval;
    ntiles::Int=16,
    min_width::Real=1e-12,
    max_halvings::Int=60,
    verbose::Bool=false,
)
    a, b = inf(Iθ), sup(Iθ)

    # θc must be strictly inside Iθ
    if !(a < θc < b)
        verbose && println("θc = $θc not strictly inside Iθ = $Iθ")
        return (:unknown, interval(θc), interval(0))
    end

    # largest symmetric radius allowed by Iθ
    rmax = min(θc - a, b - θc)
    r = rmax

    for _ in 1:max_halvings
        Ic = interval(θc, r; format=:midpoint)

        # 1) direct evaluation
        FI = f(Ic, Iκ, Iμ)
        if ispos(FI)
            verbose && println("f > 0 on ", Ic, " (direct)")
            return (:pos, Ic, FI)
        elseif isneg(FI)
            verbose && println("f < 0 on ", Ic, " (direct)")
            return (:neg, Ic, FI)
        end

        # 2) fixed-partition refinement using the generic evaluator
        FI_sub = eval_subdiv(f, Ic, Iκ, Iμ; ntiles=ntiles)
        if ispos(FI_sub)
            verbose && println("f > 0 on ", Ic, " (subdiv)")
            return (:pos, Ic, FI_sub)
        elseif isneg(FI_sub)
            verbose && println("f < 0 on ", Ic, " (subdiv)")
            return (:neg, Ic, FI_sub)
        end

        # stopping criterion
        if 2r <= min_width
            verbose && println("inconclusive at width≈", 2r, "  f(Ic)=", FI_sub)
            return (:unknown, Ic, FI_sub)
        end

        r *= 0.5
    end

    Ic = interval(θc, r; format=:midpoint)
    FI = f(Ic, Iκ, Iμ)
    return (:unknown, Ic, FI)
end

"""
    find_sign_change_centered(
        f, θc, J, Iκ, Iμ;
        min_radius=1e-12,
        max_halvings=60,
        verbose=false,
    )

Given:
  • f(Iθloc, Iκ, Iμ)        – interval evaluation of scalar f
  • θc                      – center point (approximate root)
  • J                       – θ-interval in which we search (typically a
                              certified monotonicity interval)
  • Iκ, Iμ                  – parameter intervals

Search for a *symmetric* interval [θc - r, θc + r] ⊆ J such that
f(θc - r, …) and f(θc + r, …) have **opposite strict signs**.

Algorithm:
  1) Start with the largest radius r that fits in J.
  2) Evaluate f at θL = θc - r and θR = θc + r (thin intervals).
  3) If signs are opposite → success.
  4) Otherwise halve r and repeat until r ≤ min_radius.

Returns:
    (status::Symbol, Iexist::Interval, fL::Interval, fR::Interval)

  • status = :ok       → sign change certified on Iexist = [θL, θR]
  • status = :unknown  → no sign change certified down to radius≈min_radius
"""
function find_sign_change_centered(
    f::Function,
    θc,
    J::Interval,
    Iκ::Interval,
    Iμ::Interval;
    min_radius::Real=1e-12,
    max_halvings::Int=60,
    verbose::Bool=false,
)
    a, b = inf(J), sup(J)

    # θc must lie strictly inside J
    if !(a < θc < b)
        verbose && println("θc = $θc not strictly inside J = $J")
        return (:unknown, interval(θc), interval(0), interval(0))
    end

    # largest symmetric radius allowed by J
    rmax = min(θc - a, b - θc)
    r = rmax

    # keep track of last interval / values in case we fail
    last_Iexist = interval(θc)
    last_fL = interval(0)
    last_fR = interval(0)

    for k in 1:max_halvings
        θL = θc - r
        θR = θc + r
        IL = interval(θL)
        IR = interval(θR)

        fL = f(IL, Iκ, Iμ)
        fR = f(IR, Iκ, Iμ)

        last_Iexist = interval(θL, θR)
        last_fL, last_fR = fL, fR

        left_neg_right_pos = isneg(fL) && ispos(fR)
        left_pos_right_neg = ispos(fL) && isneg(fR)

        if left_neg_right_pos || left_pos_right_neg
            verbose && println("Sign change certified at radius r = $r")
            return (:ok, last_Iexist, fL, fR)
        end

        verbose && println("No sign change at r = $r, fL=$fL, fR=$fR")

        if r <= min_radius
            verbose && println("Stopping: r ≤ min_radius ($min_radius)")
            return (:unknown, last_Iexist, last_fL, last_fR)
        end

        r *= 0.5
    end

    # Reached max_halvings without success
    verbose && println("Stopping: reached max_halvings = $max_halvings")
    return (:unknown, last_Iexist, last_fL, last_fR)
end

"""
    find_sign_change_centered_refine(
        f, θc, J, Iκ, Iμ;
        min_radius=1e-12,
        max_halvings=60,
        verbose=false,
    )

Find a symmetric sign-change bracket around θc inside J, then keep halving
to refine it. Once a sign-change has been found, STOP at the first failure
and return the last successful (smallest) bracket.

Returns:
    (status::Symbol, Iexist::Interval, fL::Interval, fR::Interval)

status:
  • :ok       → found (and refined) a sign-change bracket
  • :unknown  → never found a sign-change bracket
"""
function find_sign_change_centered_refine(
    f::Function,
    θc,
    J::Interval,
    Iκ::Interval,
    Iμ::Interval;
    min_radius::Real=1e-12,
    max_halvings::Int=100,
    verbose::Bool=false,
)
    a, b = inf(J), sup(J)

    # θc must lie strictly inside J
    if !(a < θc < b)
        verbose && println("θc = $θc not strictly inside J = $J")
        return (:unknown, interval(θc), interval(0), interval(0))
    end

    # largest symmetric radius allowed by J
    r = min(θc - a, b - θc)

    found = false
    best_Iexist = interval(θc)
    best_fL = interval(0)
    best_fR = interval(0)
    best_r = r

    for k in 1:max_halvings
        θL = θc - r
        θR = θc + r

        fL = f(interval(θL), Iκ, Iμ)
        fR = f(interval(θR), Iκ, Iμ)

        ok_sc = (isneg(fL) && ispos(fR)) || (ispos(fL) && isneg(fR))

        if ok_sc
            found = true
            best_Iexist = interval(θL, θR)
            best_fL, best_fR = fL, fR
            best_r = r
            verbose && println("✓ sign change at r = $r")
        else
            verbose && println("✗ no sign change at r = $r  fL=$fL  fR=$fR")
            if found
                verbose && println("Stopping: sign change was lost after success; returning last good r = $best_r")
                break
            end
        end

        if r <= min_radius
            break
        end
        r *= 0.5
    end

    return found ? (:ok, best_Iexist, best_fL, best_fR) :
           (:unknown, interval(θc), interval(0), interval(0))
end


"""
    validate_root_sign_mono(
        f, fθ,
        θc,
        Iθ, Iκ, Iμ;
        ntiles=16,
        min_width=1e-12,
        min_radius=1e-12,
        verbose=false,
    )

General root validation for functions f(θ, κ, μ) using:

  1. Monotonicity via fθ on a centered interval J around θc.
  2. A refined sign-change interval Iexist ⊆ J for f.

Inputs
------
  • f        : interval map  f(Iθloc, Iκ, Iμ)
  • fθ       : interval map  ∂f/∂θ(Iθloc, Iκ, Iμ)
  • θc       : approximate root (Float64 / BigFloat)
  • Iθ       : θ-interval in which we want to validate (search domain)
  • Iκ, Iμ   : parameter intervals (already built)

Keywords
--------
  • ntiles      : tiles for the subdivision in the monotonicity check
  • min_width   : minimal width for the monotonicity interval around θc
  • min_radius  : minimal radius for the sign-change search
  • verbose     : print diagnostics

Returns
-------
(ok::Bool, data::NamedTuple)

If ok == true:
    data.J       : monotonicity interval around θc
    data.Iexist  : refined existence interval with sign change
    data.fL      : f at left endpoint of Iexist
    data.fR      : f at right endpoint of Iexist
    data.Iθ      : the θ-interval used
    data.Iκ      : Iκ
    data.Iμ      : Iμ

If ok == false:
    data contains partial info and a `reason` symbol.
"""
function validate_root_sign_mono(
    f::Function,
    fθ::Function,
    θc::Real,
    Iθ::Interval,
    Iκ::Interval,
    Iμ::Interval;
    ntiles::Int=16,
    min_width::Real=1e-12,
    min_radius::Real=1e-12,
    verbose::Bool=false,
)
    a, b = inf(Iθ), sup(Iθ)

    # guard: θc must lie inside Iθ
    if !(a ≤ θc ≤ b)
        verbose && println("θc = $θc not inside Iθ = $Iθ")
        data = (Iκ=Iκ, Iμ=Iμ, Iθ=Iθ, reason=:root_outside_theta_interval)
        return false, data
    end

    # 1) monotonicity interval J via derivative
    status_J, J, _dI = certify_sign_interval_centered(
        fθ, θc, Iθ, Iκ, Iμ;
        ntiles=ntiles,
        min_width=min_width,
        max_halvings=60,
        verbose=verbose,
    )

    if status_J === :unknown
        verbose && println("validate_root_sign_mono: could not certify monotonicity near θc.")
        data = (Iκ=Iκ, Iμ=Iμ, Iθ=Iθ, J=J, reason=:no_monotonicity)
        return false, data
    end

    # 2) refined existence interval Iexist via sign change inside J
    status_sc, Iexist, fL, fR = find_sign_change_centered_refine(
        f, θc, J, Iκ, Iμ;
        min_radius=min_radius,
        max_halvings=100,
        verbose=verbose,
    )

    if status_sc !== :ok
        verbose && println("validate_root_sign_mono: no certified sign change in monotone interval.")
        data = (Iκ=Iκ, Iμ=Iμ, Iθ=Iθ, J=J, Iexist=Iexist, fL=fL, fR=fR,
            reason=:no_sign_change)
        return false, data
    end

    data = (
        Iκ=Iκ,
        Iμ=Iμ,
        Iθ=Iθ,
        J=J,
        Iexist=Iexist,
        fL=fL,
        fR=fR,
    )
    return true, data
end


##################################
### Stability validation functions
##################################

# [0,0] with the same numeric type as I
@inline zero_interval_like(I::Interval) = Interval(zero(inf(I)))

@inline function sign_tag(I::Interval)
    a, b = inf(I), sup(I)

    # guard for NaN endpoints (or anything non-finite)
    if !(isfinite(a) && isfinite(b))
        return :unknown
    end

    if a > 0
        return :pos
    elseif b < 0
        return :neg
    elseif a == 0 && b == 0
        return :zero
    elseif a ≤ 0 ≤ b
        return :mixed
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

    verbose && println("det(θ) ∈ ", Idet, " ⇒ ", tgdet)

    if tgdet === :neg
        return true, :saddle_hyperbolic, Idet,
        IntervalArithmetic.interval(0),
        IntervalArithmetic.interval(0)
    elseif tgdet !== :pos
        return false, :unknown, Idet,
        IntervalArithmetic.interval(0),
        IntervalArithmetic.interval(0)
    end

    Ilam2 = λ₂(Iθ, Iκ, Iμ)
    tgl2 = sign_tag(Ilam2)

    verbose && println("λ₂(θ) ∈ ", Ilam2, " ⇒ ", tgl2)

    if tgl2 === :pos
        return true, :lyapunov_stable, Idet, Ilam2,
        IntervalArithmetic.interval(0)
    elseif tgl2 !== :neg
        return false, :unknown, Idet, Ilam2,
        IntervalArithmetic.interval(0)
    end

    IΔ = Δ(Iθ, Iκ, Iμ)
    tgΔ = sign_tag(IΔ)

    verbose && println("Δ(θ) ∈ ", IΔ, " ⇒ ", tgΔ)

    if tgΔ === :pos
        return true, :elliptic, Idet, Ilam2, IΔ
    elseif tgΔ === :neg
        return true, :linearly_unstable, Idet, Ilam2, IΔ
    else
        return false, :unknown, Idet, Ilam2, IΔ
    end
end

@inline nonzero(I::Interval) = !(inf(I) ≤ 0 ≤ sup(I))
"""
    complements_nonzero(f, I, J, Iκ, Iμ; ntiles=16)

Return `true` iff `f(θ,κ,μ) ≠ 0` is certified on both complements
`I \\ J`.

Uses:
  • direct interval evaluation first
  • fixed subdivision fallback if needed
"""
# function complements_nonzero(
#     f::Function,
#     I::Interval,
#     J::Interval,
#     Iκ::Interval,
#     Iμ::Interval;
#     ntiles::Int=16,
# )

#     a, b = inf(I), sup(I)
#     c, d = inf(J), sup(J)

#     # require J ⊆ I
#     a ≤ c ≤ d ≤ b || return false

#     # left complement
#     if a < c
#         L = interval(a, c)
#         FL = f(L, Iκ, Iμ)
#         if !nonzero(FL)
#             FL = eval_subdiv(f, L, Iκ, Iμ; ntiles=ntiles)
#             nonzero(FL) || return false
#         end
#     end

#     # right complement
#     if d < b
#         R = interval(d, b)
#         FR = f(R, Iκ, Iμ)
#         if !nonzero(FR)
#             FR = eval_subdiv(f, R, Iκ, Iμ; ntiles=ntiles)
#             nonzero(FR) || return false
#         end
#     end

#     return true
# end


# @inline nonzero(I::Interval) = !(0 ∈ I)

"""
    complements_nonzero(f, I, J, Iκ, Iμ; ntiles=16)

Return true iff f(θ,κ,μ) ≠ 0 is certified on both complements I \\ J.

Robust to tiny rounding where J is not exactly a subset of I:
we clip J to I instead of returning false.
"""
function complements_nonzero(
    f::Function,
    I::Interval,
    J::Interval, 
    Iκ::Interval,
    Iμ::Interval;
    ntiles::Int=16,
)
    a, b = inf(I), sup(I)

    # clip J into I (instead of failing if J⊄I due to rounding)
    Jc = intersect_interval(J, I)
    isempty_interval(Jc) && return nonzero(eval_subdiv(f, I, Iκ, Iμ; ntiles=ntiles))

    c, d = inf(Jc), sup(Jc)

    # left complement [a,c]
    if a < c
        L = interval(a, c)
        FL = f(L, Iκ, Iμ)
        if !nonzero(FL)
            FL = eval_subdiv(f, L, Iκ, Iμ; ntiles=ntiles)
            nonzero(FL) || return false
        end
    end

    # right complement [d,b]
    if d < b
        R = interval(d, b)
        FR = f(R, Iκ, Iμ)
        if !nonzero(FR)
            FR = eval_subdiv(f, R, Iκ, Iμ; ntiles=ntiles)
            nonzero(FR) || return false
        end
    end

    return true
end

"""
    finalize_root_label(comp_ok; one_label=:one_root, unknown_label=:unknown)

Given:
  • comp_ok :: Bool   (output of validate_complements / complements_nonzero)

Returns:
  • (ok::Bool, label::Symbol)

Rules:
  • comp_ok == true  → one root certified
  • comp_ok == false → unknown
"""
function finalize_root_label(
    comp_ok::Bool;
    one_label::Symbol=:one_root,
    unknown_label::Symbol=:unknown,
)
    if comp_ok
        return true, one_label
    else
        return false, unknown_label
    end
end

"""
    complements_nonzero2(f, I, J1, J2, Iκ, Iμ; ntiles=16)

Return `true` iff `f(θ,κ,μ) ≠ 0` is certified on the complement
`I \\ (J1 ∪ J2)`.

Assumes `f` is interval-evaluable in θ,κ,μ.
Uses direct interval evaluation first, then subdivision fallback if needed.

Notes:
- J1 and J2 may overlap or be in any order.
- Requires both J1 ⊆ I and J2 ⊆ I (otherwise returns false).
"""
function complements_nonzero2(
    f::Function,
    I::Interval,
    J1::Interval,
    J2::Interval,
    Iκ::Interval,
    Iμ::Interval;
    ntiles::Int=16,
)
    a, b = inf(I), sup(I)

    # ensure J1,J2 are inside I
    a ≤ inf(J1) ≤ sup(J1) ≤ b || return false
    a ≤ inf(J2) ≤ sup(J2) ≤ b || return false

    # sort them by left endpoint
    if inf(J2) < inf(J1)
        J1, J2 = J2, J1
    end

    c1, d1 = inf(J1), sup(J1)
    c2, d2 = inf(J2), sup(J2)

    # merge if overlapping / touching (treat as one excluded block)
    if c2 ≤ d1
        c2 = c1
        d2 = max(d1, d2)
        c1 = c2
        d1 = d2
        # now "excluded" is [c1,d1] only
        # check I \ [c1,d1]
        # left piece
        if a < c1
            L = interval(a, c1)
            FL = f(L, Iκ, Iμ)
            if !nonzero(FL)
                FL = eval_subdiv(f, L, Iκ, Iμ; ntiles=ntiles)
                nonzero(FL) || return false
            end
        end
        # right piece
        if d1 < b
            R = interval(d1, b)
            FR = f(R, Iκ, Iμ)
            if !nonzero(FR)
                FR = eval_subdiv(f, R, Iκ, Iμ; ntiles=ntiles)
                nonzero(FR) || return false
            end
        end
        return true
    end

    # disjoint excluded intervals: J1=[c1,d1], J2=[c2,d2] with d1 < c2

    # left complement: [a, c1]
    if a < c1
        L = interval(a, c1)
        FL = f(L, Iκ, Iμ)
        if !nonzero(FL)
            FL = eval_subdiv(f, L, Iκ, Iμ; ntiles=ntiles)
            nonzero(FL) || return false
        end
    end

    # middle complement: [d1, c2]
    if d1 < c2
        M = interval(d1, c2)
        FM = f(M, Iκ, Iμ)
        if !nonzero(FM)
            FM = eval_subdiv(f, M, Iκ, Iμ; ntiles=ntiles)
            nonzero(FM) || return false
        end
    end

    # right complement: [d2, b]
    if d2 < b
        R = interval(d2, b)
        FR = f(R, Iκ, Iμ)
        if !nonzero(FR)
            FR = eval_subdiv(f, R, Iκ, Iμ; ntiles=ntiles)
            nonzero(FR) || return false
        end
    end

    return true
end



"""
    classify_center_equilibrium(Ix0, Iκ, Iμ; verbose=false)

Interval implementation of the center-equilibrium test:

Define
    Λ  = ( sin(√κ q₁) / ( μ sin(√κ q₂) ) )^(1/3)
    d₁ = asin( Λ * sin(x0) )
    d₂ = x0

Steps:
1) If det(d₁,d₂) < 0  → :center_saddle_unstable
2) else if tr(d₁,d₂) > 0 → :lyapunov_stable
3) else if C(d₁,d₂) > 0 → :elliptic
   else if C(d₁,d₂) < 0 → :complex_saddle_unstable
   else → :unknown

Returns:
    (ok::Bool, label::Symbol, Idet, Itr, IC, Id1, Id2, IΛ)
"""
function classify_center_equilibrium(Ix0::Interval, Iκ::Interval, Iμ::Interval; verbose=false)

    oneI = interval(1)
    twoI = interval(2)
    threeI = interval(3)
    fourI = interval(4)
    sixI = interval(6)

    r = sqrt(Iκ)
    q1 = q₁(Iκ, Iμ)
    q2 = q₂(Iκ, Iμ)

    # Λ = ( sin(√κ q1) / ( μ sin(√κ q2) ) )^(1/3)
    IΛ = cbrt(sin(r * q1) / (Iμ * sin(r * q2)))

    # d1 = asin( Λ * sin(x0) ), d2 = x0
    arg = IΛ * sin(Ix0)
    if !(inf(arg) ≥ -1 && sup(arg) ≤ 1)
        verbose && println("asin argument not in [-1,1]: ", arg)
        z = interval(0)
        return false, :unknown, z, z, z, z, Ix0, IΛ
    end

    Id1 = asin(arg)
    Id2 = Ix0

    # ---- det(d1,d2) ----------------------------------------------------

    sr = sin(r)
    sr4 = sr^4

    s_rq1 = sin(r * q1)
    s_rq2 = sin(r * q2)

    pref = (threeI * s_rq1 * s_rq2) / sr4

    B1 = s_rq2 * cos(Id1) + s_rq1 * cos(Id2)

    B2 = (oneI + twoI * cos(twoI * Id2)) * s_rq1 * cos(Id1) +
         (oneI + twoI * cos(twoI * Id1)) * s_rq2 * cos(Id2)

    Idet = pref * B1 * B2

    # ---- tr(d1,d2) -----------------------------------------------------

    Itr = (oneI / (sr^2)) * (
        (oneI - fourI * cos(Id2)^2) * s_rq1^2
        -
        sixI * cos(Id1) * cos(Id2) * s_rq1 * s_rq2
        +
        (oneI - fourI * cos(Id1)^2) * s_rq2^2
    )

    # ---- C(d1,d2) ------------------------------------------------------

    Im = -(sin(Id1)^2 * sin(Id2)^2) /
         ((cos(r) - cos(Id1 + Id2)) * (cos(r) - cos(Id1 - Id2)))

    IC = Im - fourI * Idet

    # ---- Classification ------------------------------------------------

    tgdet = sign_tag(Idet)
    verbose && println("det ∈ ", Idet, " ⇒ ", tgdet)

    if tgdet === :neg
        return true, :center_saddle_unstable, Idet, Itr, IC, Id1, Id2, IΛ
    elseif tgdet !== :pos
        return false, :unknown, Idet, Itr, IC, Id1, Id2, IΛ
    end

    tgtr = sign_tag(Itr)
    verbose && println("tr ∈ ", Itr, " ⇒ ", tgtr)

    if tgtr === :pos
        return true, :lyapunov_stable, Idet, Itr, IC, Id1, Id2, IΛ
    elseif tgtr !== :neg
        return false, :unknown, Idet, Itr, IC, Id1, Id2, IΛ
    end

    tgC = sign_tag(IC)
    verbose && println("C ∈ ", IC, " ⇒ ", tgC)

    if tgC === :pos
        return true, :elliptic, Idet, Itr, IC, Id1, Id2, IΛ
    elseif tgC === :neg
        return true, :complex_saddle_unstable, Idet, Itr, IC, Id1, Id2, IΛ
    else
        return false, :unknown, Idet, Itr, IC, Id1, Id2, IΛ
    end
end