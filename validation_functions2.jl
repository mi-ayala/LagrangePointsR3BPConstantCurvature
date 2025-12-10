using IntervalArithmetic

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

using IntervalArithmetic

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

using IntervalArithmetic

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
    status_sc, Iexist, fL, fR = find_sign_change_centered(
        f, θc, J, Iκ, Iμ;
        min_radius=min_radius,
        max_halvings=60,
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