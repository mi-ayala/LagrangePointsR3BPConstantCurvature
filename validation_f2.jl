using IntervalArithmetic
using Printf
using Roots
using CairoMakie

include("functions_helpers.jl")
include("functions_stability.jl")
include("functions_validation.jl")
include("functions_f2.jl")
include("functions_plots.jl")


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

    Ileft = sqrt(Iκ) * q₂(Iκ, Iμ)          # interval lower bound
    Iright = interval(pi / 2)                 # thin interval upper bound
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
    tg = sign_tag(If2)

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
    rκ = Δκ / 2
    rμ = Δμ / 2

    boxes = Tuple{Symbol,Float64,Float64,Float64,Float64,Bool,Bool,Symbol,Symbol}[]
    total = nx * ny
    seeded = 0
    validated = 0

    # iterate cell centers
    for i in 1:nx, j in 1:ny
        κc = κ_lo + (i - 0.5) * Δκ
        μc = μ_lo + (j - 0.5) * Δμ

        (θmin, θmax) = Intθ_2(κc, μc)

        br = brackets_signchanges(θ -> f₂θ(θ, κc, μc), θmin, θmax)

        if isempty(br)
            continue
        end
        seeded += 1 # Seeded means we found a root 

        θM = θM_from_br(br, κc, μc; atol=1e-12)
        # 2) refine the root 
        root = find_zero(θ -> f₂θ(θ, κc, μc), θM; atol=1e-12)

        Iκ = interval(κc, rκ; format=:midpoint)
        Iμ = interval(μc, rμ; format=:midpoint)

        θmax_refined = sup(Intθ_2(Iκ, Iμ)[2])
        θmin_refined = inf(Intθ_2(Iκ, Iμ)[1])

        ok_root, data = validate_root_sign_mono(
            f₂θ,
            f₂θθ,
            root,
            interval(θmin_refined, θmax_refined),
            Iκ,
            Iμ;
            ntiles=16,
            min_width=1e-12,
            min_radius=1e-12,
            verbose=false,
        )

        keep_invalid = true

        if ok_root

            clasif = classify_f2_rootcount(data.Iexist, Iκ, Iμ; verbose=false)

            if clasif[2] === :unknown
                push!(boxes, (:unknown, κc, μc, rκ, rμ, false, false, :unknown, :unknown))

            elseif clasif[2] === :zero_roots
                push!(boxes, (:zero_roots, κc, μc, rκ, rμ, false, false, :unknown, :unknown))
                validated += 1

            elseif clasif[2] === :two_roots

                br = brackets_signchanges(θ -> f₂(θ, κc, μc), θmin, θmax)

                if length(br) < 2
                    push!(boxes, (:two_roots, κc, μc, rκ, rμ, false, false, :unknown, :unknown))
                    continue
                end

                root1 = find_zero(θ -> f₂(θ, κc, μc), br[1]; atol=1e-12)
                root2 = find_zero(θ -> f₂(θ, κc, μc), br[2]; atol=1e-12)
                
                ### order roots by absolute value
                if abs(root2) < abs(root1)
                    (root1, root2) = (root2, root1)
                end

                ok_root1, data1 = validate_root_sign_mono(
                    f₂,
                    f₂θ,
                    root1,
                    interval(θmin_refined, θmax_refined),
                    Iκ,
                    Iμ;
                    ntiles=30,
                    min_width=1e-12,
                    min_radius=1e-12,
                    verbose=false,
                )

                ok_root2, data2 = validate_root_sign_mono(
                    f₂,
                    f₂θ,
                    root2,
                    interval(θmin_refined, θmax_refined),
                    Iκ,
                    Iμ;
                    ntiles=30,
                    min_width=1e-12,
                    min_radius=1e-12,
                    verbose=false,
                )

                _, stability_label1 = classify_equilibrium(data1.Iexist, Iκ, Iμ)
                _, stability_label2 = classify_equilibrium(data2.Iexist, Iκ, Iμ)




                push!(boxes, (clasif[2], κc, μc, rκ, rμ, ok_root1, ok_root2, stability_label1, stability_label2))
                validated += 1
            end


        elseif keep_invalid

            push!(boxes, (:unknown, κc, μc, rκ, rμ, false, false, :unknown, :unknown))
        end

    end

    if verbose
        println("cover_rectangle_grid_rect:")
        println("  Δκ=$(Δκ), Δμ=$(Δμ); rκ=$(rκ), rμ=$(rμ)")
        println("  cells total=$(total), seeded=$(seeded), validated=$(validated)")
    end

    stats = (Δκ=Δκ, Δμ=Δμ, rκ=rκ, rμ=rμ, total=total, seeded=seeded, validated=validated)
    return boxes, stats
end


κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 3 * 2464, 3 * 89
boxes, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


# --- Root count plot ---
fig1, ax1 = plot_regions_all_labels_from_boxes_quick(boxes;
    label_type = :root,
    title = "Root-count",
    alpha = 0.9,
)

save("figures/f2_root_regions_quick.png", fig1)

# --- Stability plot (root 1) ---
fig2, ax2 = plot_regions_all_labels_from_boxes_quick(boxes;
    label_type = :stability,
    root_index = 1,
    alpha = 0.9,
    title = "Stability (root 1)",
)

save("figures/f2_stability_root1_quick.png", fig2)

# --- Stability plot (root 2) ---
fig3, ax3 = plot_regions_all_labels_from_boxes_quick(boxes;
    label_type = :stability,
    root_index = 2,
    alpha = 0.9,
    title = "Stability (root 2)",
)

save("figures/f2_stability_root2_quick.png", fig3)




# # --- Root count plot ---
# fig1, ax1 = plot_regions_all_labels_from_boxes(boxes;
#     label_type = :root,
#     title = "Root-count",
# )

# save("figures/f2_root_regions_full.png", fig1)

# # --- Stability plot (root 1) ---
# fig2, ax2 = plot_regions_all_labels_from_boxes(boxes;
#     label_type = :stability,
#     root_index = 1,
#     title = "Stability (root 1)",
# )

# save("figures/f2_stability_root1.png", fig2)

# # --- Stability plot (root 2) ---
# fig3, ax3 = plot_regions_all_labels_from_boxes(boxes;
#     label_type = :stability,
#     root_index = 2,
#     title = "Stability (root 2)",
# )

# save("figures/f2_stability_root2.png", fig3)

