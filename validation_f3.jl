using IntervalArithmetic
using Printf
using Roots
using CairoMakie

include("functions_helpers.jl")
include("functions_validation.jl")
include("functions_f3.jl")
include("functions_plots.jl")


function classify_f3_rootcount(Iθm::Interval, Iκ::Interval, Iμ::Interval; verbose::Bool=false)

    Ileft = interval(-pi / 2)
    Iright = -(sqrt(Iκ) * q₁(Iκ, Iμ))
    Iallowed = interval(inf(Ileft), sup(Iright))

    inside = (inf(Iθm) ≥ inf(Iallowed)) && (sup(Iθm) ≤ sup(Iallowed))

    if verbose
        println("θ_m ∈ ", Iθm, " ; allowed ∈ ", Iallowed, " ; inside? ", inside)
    end

    if !inside

        return true, :zero_roots, Iθm, interval(0), Iallowed
    end


    If3 = f₃(Iθm, Iκ, Iμ)
    tg = sign_tag(If3)

    if verbose
        println("f₃(θ_m) ∈ ", If3, " ⇒ ", tg)
    end

    if tg === :neg
        return true, :two_roots, Iθm, If3, Iallowed
    elseif tg === :pos
        return true, :zero_roots, Iθm, If3, Iallowed
    else
        return false, :unknown, Iθm, If3, Iallowed
    end
end

function θm_from_br(br, κ, μ; atol=1e-12, verbose=false)
    bestθ = NaN
    bestval = Inf
    found = false

    for x0 in br
        θ = find_zero(θ -> f₃θ(θ, κ, μ), x0; atol=atol)

        s2 = f₃θθ(θ, κ, μ)      # second derivative test
        if s2 > 0               # local minimum
            val = f₃(θ, κ, μ)   # height of the minimum
            if !found || val < bestval
                bestval = val
                bestθ = θ
                found = true
            end
        end

        verbose && println("x0=$x0  θ*=$θ  f₃θθ=$s2  f₃=$(f₃(θ, κ, μ))")
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

        (θmin, θmax) = Intθ_3(κc, μc)

        br = brackets_signchanges(θ -> f₃θ(θ, κc, μc), θmin, θmax)

        if isempty(br)
            continue
        end
        seeded += 1 # Seeded means we found a root 

        θm = θm_from_br(br, κc, μc; atol=1e-12)
        # 2) refine the root 
        root = find_zero(θ -> f₃θ(θ, κc, μc), θm; atol=1e-12)

        Iκ = interval(κc, rκ; format=:midpoint)
        Iμ = interval(μc, rμ; format=:midpoint)

        θmax_refined = sup(Intθ_3(Iκ, Iμ)[2])
        θmin_refined = inf(Intθ_3(Iκ, Iμ)[1])

        ok_root, data = validate_root_sign_mono(
            f₃θ,
            f₃θθ,
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

            clasif = classify_f3_rootcount(data.Iexist, Iκ, Iμ; verbose=false)

            if clasif[2] === :unknown
                push!(boxes, (:unknown, κc, μc, rκ, rμ, false, false, :unknown, :unknown))

            elseif clasif[2] === :zero_roots
                push!(boxes, (:zero_roots, κc, μc, rκ, rμ, false, false, :unknown, :unknown))
                validated += 1  

            else

                br = brackets_signchanges(θ -> f₃(θ, κc, μc), θmin, θmax)


                root1 = find_zero(θ -> f₃(θ, κc, μc), br[1]; atol=1e-12)
                root2 = find_zero(θ -> f₃(θ, κc, μc), br[2]; atol=1e-12)

                ### order roots by absolute value
                if abs(root2) < abs(root1)
                    (root1, root2) = (root2, root1)
                end

                ok_root1, data1 = validate_root_sign_mono(
                    f₃,
                    f₃θ,
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
                    f₃,
                    f₃θ,
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




                push!(boxes, ( clasif[2], κc, μc, rκ, rμ, ok_root1, ok_root2, stability_label1, stability_label2))
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

nx, ny = 2 * 2464, 2 * 89
boxes, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)





# --- Root count plot ---
fig1, ax1 = plot_regions_all_labels_from_boxes_quick(boxes;
    label_type = :root,
    title = "Root-count",
)

save("figures/f3_root_regions_quick.png", fig1)

# --- Stability plot (root 1) ---
fig2, ax2 = plot_regions_all_labels_from_boxes_quick(boxes;
    label_type = :stability,
    root_index = 1,
    title = "Stability (root 1)",
)

save("figures/f3_stability_root1_quick.png", fig2)

# --- Stability plot (root 2) ---
fig3, ax3 = plot_regions_all_labels_from_boxes_quick(boxes;
    label_type = :stability,
    root_index = 2,
    title = "Stability (root 2)",
)

save("figures/f3_stability_root2_quick.png", fig3)





# fig1, ax = plot_regions_all_labels_from_boxes(
#     boxes;
#     label_type=:root,
#     κlims=(0.0, (0.5 * pi)^2),
#     μlims=(0.0, 0.0886765),
#     resolution=2000,
#     alpha=0.7,
#     show_boundary=true,
#     boundary_lw=2.0,
#     title="Root-count Regions",
# )


# CairoMakie.save("figures/f3_root_regions.png", fig1)
