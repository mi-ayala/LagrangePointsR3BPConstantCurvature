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

        (θmin, θmax) = Intθ_1(κc, μc)

        br = brackets_signchanges(θ -> f₁(θ, κc, μc), θmin, θmax)

        if isempty(br)
            continue
        end
        seeded += 1

        root = find_zero(θ -> f₁(θ, κc, μc), br[1]; atol=1e-12)

        Iκ = interval(κc, rκ; format=:midpoint)
        Iμ = interval(μc, rμ; format=:midpoint)

        θmax_refined = sup(θ1_max(Iκ, Iμ))

        ok_root, data = validate_root_sign_mono(
            f₁,
            f₁θ,
            root,
            interval(0.0, θmax_refined),
            Iκ,
            Iμ;
            ntiles=16,
            min_width=1e-12,
            min_radius=1e-12,
            verbose=false,
        )

        keep_invalid = true


        ### Domain checking
        (I_left,I_right) = Intθ_1(Iκ, Iμ)
        if sup(I_left) > inf(data.Iθ) || inf(I_right) < sup(data.Iθ)
            ok_root = false
        end


        if ok_root

            _, root_label = finalize_root_label(complements_nonzero(f₁, interval(0.0, θmax_refined), data.J, Iκ, Iμ))

            _, stability_label = classify_equilibrium(data.Iexist, Iκ, Iμ)

            if root_label === :unknown

                push!(boxes, (root_label, κc, μc, rκ, rμ, ok_root, false, :unknown, :zero_roots))

            else
                push!(boxes, (root_label, κc, μc, rκ, rμ, ok_root, false, stability_label, :zero_roots))
                validated += 1
            end


        elseif keep_invalid
            push!(boxes, (:unknown, κc, μc, rκ, rμ, ok_root, false, :unknown, :zero_roots))
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
μ_lo, μ_hi = 0.02, 0.0886765

nx, ny = 4*2464,  4*89
boxes1_validated, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)

κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.02

nx, ny = 4* 2464, 4 * 89
boxes1_validated2, stats2 = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)

append!(boxes1_validated, boxes1_validated2)

fig1, ax1 = plot_regions_all_labels_from_boxes(
    boxes1_validated;
    label_type=:root,
    overlay_rect=(0.0, pi^2 / 16, 0.0, 0.0886765),
    overlay_color=RGBAf(0.97, 0.98, 1.00, 1.0)
)
save("figures/validated/f1_root_regions.png", fig1)

fig2, ax2 = plot_regions_all_labels_from_boxes(boxes1_validated;
    label_type=:stability,
    root_index=1,
    overlay_rect=(0.0, pi^2 / 16, 0.0, 0.0886765),
    overlay_color=RGBAf(0.97, 0.98, 1.00, 1.0),
)
save("figures/validated/f1_stability.png", fig2)

