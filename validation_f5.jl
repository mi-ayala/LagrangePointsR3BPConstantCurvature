using Roots
using CairoMakie
using IntervalArithmetic
using RC3BP


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

        (θmin, θmax) = Intθ_5(κc, μc)

        br = brackets_signchanges(θ -> f₅(θ, κc, μc), θmin, θmax)

        Iκ = interval(κc, rκ; format=:midpoint)
        Iμ = interval(μc, rμ; format=:midpoint)

        θmax_refined = sup(Intθ_5(Iκ, Iμ)[2])
        θmin_refined = inf(Intθ_5(Iκ, Iμ)[1])

        # ------------------------
        # Case 0: no brackets
        # ------------------------
        if isempty(br)
            # certify “no root” by proving f5 has fixed sign on the whole θ-interval
            FI = eval_subdiv(f₅, interval(θmin_refined, θmax_refined), Iκ, Iμ; ntiles=50)

            if ispos(FI) || isneg(FI)
                push!(boxes, (:zero_roots, κc, μc, rκ, rμ, false, false, :zero_roots, :zero_roots))
                validated += 1
            else
                push!(boxes, (:unknown, κc, μc, rκ, rμ, false, false, :unknown, :unknown))
            end

        # ------------------------
        # Case 1: exactly one bracket
        # ------------------------
        elseif length(br) == 1
            root = find_zero(θ -> f₅(θ, κc, μc), br[1]; atol=1e-12)

            ok_root, data = validate_root_sign_mono(
                f₅,
                f₅θ,
                root,
                interval(θmin_refined, θmax_refined),
                Iκ,
                Iμ;
                ntiles=60,
                min_width=1e-12,
                min_radius=1e-12,
                verbose=false,
            )

            if ok_root
                comp_ok = complements_nonzero(
                    f₅,
                    interval(θmin_refined, θmax_refined),
                    data.J,
                    Iκ,
                    Iμ,
                )

                if comp_ok
                    _, stability_label = classify_center_equilibrium(data.Iexist, Iκ, Iμ)
                    push!(boxes, (:one_root, κc, μc, rκ, rμ, true, false, stability_label, :one_root))
                    validated += 1
                else
                    push!(boxes, (:unknown, κc, μc, rκ, rμ, true, false, :unknown, :unknown))
                end
            else
                push!(boxes, (:unknown, κc, μc, rκ, rμ, false, false, :unknown, :unknown))
            end

        # ------------------------
        # Case 2: two or more brackets (attempt “two roots”)
        # ------------------------
        else
            root1 = find_zero(θ -> f₅(θ, κc, μc), br[1]; atol=1e-12)
            root2 = find_zero(θ -> f₅(θ, κc, μc), br[2]; atol=1e-12)

            # order roots by absolute value
            if abs(root2) < abs(root1)
                (root1, root2) = (root2, root1)
            end

            ok_root1, data1 = validate_root_sign_mono(
                f₅,
                f₅θ,
                root1,
                interval(θmin_refined, θmax_refined),
                Iκ,
                Iμ;
                ntiles=40,
                min_width=1e-12,
                min_radius=1e-12,
                verbose=false,
            )
            ok_root2, data2 = validate_root_sign_mono(
                f₅,
                f₅θ,
                root2,
                interval(θmin_refined, θmax_refined),
                Iκ,
                Iμ;
                ntiles=40,
                min_width=1e-12,
                min_radius=1e-12,
                verbose=false,
            )

            if ok_root1 && ok_root2
                _, stability_label1 = classify_center_equilibrium(data1.Iexist, Iκ, Iμ)
                _, stability_label2 = classify_center_equilibrium(data2.Iexist, Iκ, Iμ)

                comp_ok = complements_nonzero2(
                    f₅,
                    interval(θmin_refined, θmax_refined),
                    data1.J,
                    data2.J,
                    Iκ,
                    Iμ,
                )

                if comp_ok
                    push!(boxes, (:two_roots, κc, μc, rκ, rμ, true, true, stability_label1, stability_label2))
                    validated += 1
                else
                    push!(boxes, (:unknown, κc, μc, rκ, rμ, true, true, :unknown, :unknown))
                end
            else
                push!(boxes, (:unknown, κc, μc, rκ, rμ, ok_root1, ok_root2, :unknown, :unknown))
            end
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

nx, ny =  4*2464, 4*89
@time boxes5_validated, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)














# --- Root count plot ---
fig1, ax1 = plot_regions_all_labels_from_boxes(boxes5_validated;
    label_type=:root,
    # title="Root-count",
    # κlims=(0.0, 0.24),
)

save("/Users/miguelayala/Research/R3BP-on-spaces-with-constant-curvature/figures Miguel/validated/f5_root_regions.png", fig1)

# --- Stability plot (root 1) ---

fig2, ax2 = plot_regions_all_labels_from_boxes(boxes5_validated;
    label_type=:stability,
    root_index=1,
    # title="Stability (root 1)",
)

save("/Users/miguelayala/Research/R3BP-on-spaces-with-constant-curvature/figures Miguel/validated/f5_stability_root1.png", fig2)

# --- Stability plot (root 2) ---
fig3, ax3 = plot_regions_all_labels_from_boxes(boxes5_validated;
    label_type=:stability,
    root_index=2,
    # title="Stability (root 2)",
)

save("/Users/miguelayala/Research/R3BP-on-spaces-with-constant-curvature/figures Miguel/validated/f5_stability_root2.png", fig3)

