using IntervalArithmetic
using Plots
using PlotlyJS
using Roots
plotlyjs()

include("functions.jl")
include("validation_functions.jl")

κ = 2.0
μ = 0.05

# helpers for interval sign tagging
sign_tag(I) = ispos(I) ? :pos :
              isneg(I) ? :neg :
              (0 ∈ I ? :zero_or_mixed : :unknown)

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
    Idet  = detλ(Iθ, Iκ, Iμ)
    tgdet = sign_tag(Idet)

    if verbose
        println("det(θ) ∈ ", Idet, "  ⇒ ", tgdet)
    end

    if tgdet === :neg
        return true, :saddle_hyperbolic, Idet, IntervalArithmetic.Interval(0), IntervalArithmetic.Interval(0)
    end

    Ilam2  = λ₂(Iθ, Iκ, Iμ)
    tgl2   = sign_tag(Ilam2)
    if verbose
        println("λ₂(θ) ∈ ", Ilam2, "  ⇒ ", tgl2)
    end

    if tgl2 === :pos
        return true, :lyapunov_stable, Idet, Ilam2, IntervalArithmetic.Interval(0)
    end

    IΔ   = Δ(Iθ, Iκ, Iμ)
    tgΔ  = sign_tag(IΔ)
    if verbose
        println("Δ(θ) ∈ ", IΔ, "  ⇒ ", tgΔ)
    end

    if tgΔ === :pos
        return true, :elliptic, Idet, Ilam2, IΔ
    elseif tgΔ === :neg
        return true, :linearly_unstable, Idet, Ilam2, IΔ
    else
        return false, :unknown, Idet, Ilam2, IΔ
    end
end
