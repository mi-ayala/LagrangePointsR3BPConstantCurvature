using IntervalArithmetic
using Plots
using Roots
using Statistics


include("functions.jl")
include("stability_functions.jl")
include("functions_f2.jl")
include("validation_functions2.jl")

κ = 0.12
μ = 0.005

(θmin, θmax) = Intθ_2(κ, μ)

br = brackets_signchanges(θ -> f₂(θ, κ, μ), θmin, θmax)

root = find_zero(θ -> f₂(θ, κ, μ), br[1]; atol=1e-12)
root2 = find_zero(θ -> f₂(θ, κ, μ), br[2]; atol=1e-12)

Iκ = interval(κ, 1e-5; format=:midpoint)
Iμ = interval(μ, 1e-5; format=:midpoint)

ok, data = validate_root_sign_mono(
    f₂,
    f₂θ,
    root,
    interval(θmin, θmax),
    Iκ,
    Iμ;
    ntiles=16,
    min_width=1e-12,
    min_radius=1e-12,
    verbose=true,
)

ok, data2 = validate_root_sign_mono(
    f₂,
    f₂θ,
    root2,
    interval(θmin, θmax),
    Iκ,
    Iμ;
    ntiles=16,
    min_width=1e-12,
    min_radius=1e-12,
    verbose=true,
)

classify_equilibrium(data.Iexist, Iκ, Iμ)
classify_equilibrium(data2.Iexist, Iκ, Iμ)


"""
    classify_f2_rootcount(IθM, Iκ, Iμ; verbose=false)

Interval version of the 3-step test for f₂:

Step 2: check θ_M ∈ [ √κ q₂(κ,μ),  π/2 ].
        If not certified inside → :zero_roots.
Step 3: check sign of f₂(θ_M).
        If strictly positive → :two_roots.
        If strictly negative → :zero_roots.
        Else → :unknown.

Returns:
  (ok::Bool, label::Symbol, IθM, If2, Iallowed)

`ok == true` iff label ∈ (:two_roots, :zero_roots) is certified.
"""
function classify_f2_rootcount(IθM::Interval, Iκ::Interval, Iμ::Interval; verbose::Bool=false)
    # --- Step 2: θ_M must be inside [ √κ q₂, π/2 ] ------------------------

    Ileft  = sqrt(Iκ) * q₂(Iκ, Iμ)          # interval lower bound
    Iright = interval(pi/2)                 # thin interval upper bound
    Iallowed = interval(inf(Ileft), sup(Iright))

    inside = (inf(IθM) ≥ inf(Iallowed)) && (sup(IθM) ≤ sup(Iallowed))

    if verbose
        println("θ_M ∈ ", IθM, " ; allowed ∈ ", Iallowed, " ; inside? ", inside)
    end

    if !inside
        # certified outside ⇒ zero roots
        return true, :zero_roots, IθM, interval(0), Iallowed
    end

    # --- Step 3: sign of f₂ at θ_M ---------------------------------------

    If2 = f₂(IθM, Iκ, Iμ)
    tg  = sign_tag(If2)

    if verbose
        println("f₂(θ_M) ∈ ", If2, " ⇒ ", tg)
    end

    if tg === :pos
        return true, :two_roots, IθM, If2, Iallowed
    elseif tg === :neg
        return true, :zero_roots, IθM, If2, Iallowed
    else
        return false, :unknown, IθM, If2, Iallowed
    end
end

using Roots

"""
    θM_from_br(br, κ, μ; atol=1e-12)

Given `br` = vector of numeric guesses for roots of f₂θ(θ,κ,μ),
refine each by `find_zero(..., br[i]; atol=atol)`, keep those with
f₂θθ(θ*) < 0 (local maxima), and return the maximizer θ_M.
Returns NaN if no maximum is found.
"""
function θM_from_br(br, κ, μ; atol=1e-12, verbose=false)
    bestθ = NaN
    bestval = -Inf
    found = false

    for x0 in br
        θ = find_zero(θ -> f₂θ(θ, κ, μ), x0; atol=atol)

        s2 = f₂θθ(θ, κ, μ)     # second derivative test
        if s2 < 0              # local maximum
            val = f₂(θ, κ, μ)  # height of the maximum
            if !found || val > bestval
                bestval = val
                bestθ = θ
                found = true
            end
        end

        verbose && println("x0=$x0  θ*=$θ  f₂θθ=$s2  f₂=$((f₂(θ, κ, μ)))")
    end

    return found ? bestθ : NaN
end

br = brackets_signchanges(θ -> f₂θ(θ, κ, μ), θmin, θmax)


θM = θM_from_br(br, κ, μ; atol=1e-12)

root = find_zero(θ -> f₂θ(θ, κ, μ), br[1]; atol=1e-12)
root2 = find_zero(θ -> f₂θ(θ, κ, μ), br[2]; atol=1e-12)
root3 = find_zero(θ -> f₂θ(θ, κ, μ), br[3]; atol=1e-12)

ok, data = validate_root_sign_mono(
    f₂θ,
    f₂θθ,
    root,
    interval(θmin, θmax),
    Iκ,
    Iμ;
    ntiles=16,
    min_width=1e-12,
    min_radius=1e-12,
    verbose=true,
)



# f₂θθ(root, κ, μ)
# f₂θθ(root2, κ, μ)
# f₂θθ(root3, κ, μ)

# f₂(root, κ, μ)
# f₂(root2, κ, μ)
# f₂(root3, κ, μ)

classif = classify_f2_rootcount(data.Iexist, Iκ, Iμ; verbose=true)