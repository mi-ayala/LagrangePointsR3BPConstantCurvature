using Roots
using CairoMakie


include("functions_helpers.jl")
include("functions_stability.jl")
include("functions_validation.jl")
include("functions_f2.jl")
include("functions_plots.jl")
include("functions_numerics.jl")


"""
    classify_f2_rootcount(θM, κ, μ; verbose=false)

NUMERICAL (Float) version of the f₂ root-count test:

Step 2: check θ_M ∈ [ √κ q₂(κ,μ),  π/2 ].
        If not inside → :zero_roots.
Step 3: check sign of f₂(θ_M).
        If f₂(θ_M) > 0 → :two_roots.
        If f₂(θ_M) < 0 → :zero_roots.
        Else → :unknown.

Returns:
  (ok::Bool, label::Symbol, θM::Float64, f2val::Float64, allowed::Tuple{Float64,Float64})

`ok == true` iff label ∈ (:two_roots, :zero_roots) is decided from strict inequalities.
"""
function classify_f2_rootcount(θM::Real, κ::Real, μ::Real; verbose::Bool=false)
    θM = float(θM)
    κ = float(κ)
    μ = float(μ)

    left = sqrt(κ) * float(q₂(κ, μ))
    right = (pi / 2)
    allowed = (left, right)

    inside = (left ≤ θM ≤ right)

    if verbose
        println("θ_M = ", θM, " ; allowed = ", allowed, " ; inside? ", inside)
    end

    if !inside
        return true, :zero_roots, θM, 0.0, allowed
    end

    f2val = float(f₂(θM, κ, μ))

    if verbose
        println("f₂(θ_M) = ", f2val)
    end

    if f2val > 0
        return true, :two_roots, θM, f2val, allowed
    elseif f2val < 0
        return true, :zero_roots, θM, f2val, allowed
    else
        return false, :unknown, θM, f2val, allowed
    end
end


"""
    θM_from_br(br, κ, μ; atol=1e-12, verbose=false)

NUMERICAL (Float) version:

Given `br` = vector of numeric guesses for stationary points (roots of f₂θ),
refine each with `find_zero`, keep those with f₂θθ(θ*) < 0 (local maxima),
and return the maximizer θ_M (argmax of f₂ at those maxima).

Returns θ_M::Float64, or NaN if no maximum is found.
"""
function θM_from_br(br, κ, μ; atol=1e-12, verbose::Bool=false)
    κ = float(κ)
    μ = float(μ)

    bestθ = NaN
    bestval = -Inf
    found = false

    for x0 in br
        x0 = float(x0)

        θ = find_zero(θ -> f₂θ(θ, κ, μ), x0; atol=atol)

        s2 = float(f₂θθ(θ, κ, μ))   # second derivative test
        if s2 < 0                  # local maximum
            val = float(f₂(θ, κ, μ))
            if !found || val > bestval
                bestval = val
                bestθ = θ
                found = true
            end
        end

        verbose && println("x0=$x0  θ*=$θ  f₂θθ=$s2  f₂=$(float(f₂(θ, κ, μ)))")
    end

    return found ? bestθ : NaN
end





function cover_rectangle_grid_rect(κ_lo::Float64, κ_hi::Float64,
    μ_lo::Float64, μ_hi::Float64;
    nx::Int, ny::Int,
    verbose::Bool=true)

    @assert κ_hi > κ_lo "κ_hi must be > κ_lo"
    @assert μ_hi > μ_lo "μ_hi must be > μ_lo"
    @assert nx > 0 && ny > 0 "nx, ny must be positive"

    Lκ = κ_hi - κ_lo
    Lμ = μ_hi - μ_lo
    Δκ = Lκ / nx
    Δμ = Lμ / ny


    boxes = Tuple{Symbol,Float64,Float64,Float64,Float64,Symbol,Symbol}[]

    # iterate cell centers
    for i in 1:nx, j in 1:ny
        κc = κ_lo + (i - 0.5) * Δκ
        μc = μ_lo + (j - 0.5) * Δμ

        θmin = Intθ_2(κc, μc)[1]
        θmax = Intθ_2(κc, μc)[2]

        br = brackets_signchanges(θ -> f₂(θ, κc, μc), θmin, θmax)

        num_roots = length(br)

        if num_roots == 0
            push!(boxes, (:zero_roots, κc, μc, NaN, NaN, :zero_roots, :zero_roots))

        elseif num_roots == 1

            root = find_zero(θ -> f₂(θ, κc, μc), br[1]; atol=1e-12)
            stability_root = classify_equilibrium_f(root, κc, μc; verbose=false)
            push!(boxes, (:one_root, κc, μc, root, NaN, stability_root, :zero_roots))

        elseif num_roots == 2

            root1 = find_zero(θ -> f₂(θ, κc, μc), br[1]; atol=1e-12)
            stability_root1 = classify_equilibrium_f(root1, κc, μc; verbose=false)


            root2 = find_zero(θ -> f₂(θ, κc, μc), br[2]; atol=1e-12)
            stability_root2 = classify_equilibrium_f(root2, κc, μc; verbose=false)

            push!(boxes, (:two_roots, κc, μc, root1, root2, stability_root1, stability_root2))
        else
            push!(boxes, (:multiple_roots, κc, μc, NaN, NaN, :unknown, :unknown))


        end


    end


    return boxes
end


κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 100, 100
boxes = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)

# --- Root count plot ---
fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type = :root,
    ms = 4,
)

save("figures/f2_root_regions_numerics.png", fig1)

# --- Stability plot (root 1) ---
fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type = :stability,
    root_index = 1,
        ms = 4
)

save("figures/f2_stability_root1_numerics.png", fig2)

# --- Stability plot (root 2) ---
fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type = :stability,
    root_index = 2,
        ms = 4
)

save("figures/f2_stability_root2_numerics.png", fig3)




