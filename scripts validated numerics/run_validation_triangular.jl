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

        (θmin, θmax) = Intθ_g(κc, μc)

        br = brackets_signchanges(θ -> g(θ, κc, μc), θmin, θmax)

        Iκ = interval(κc, rκ; format=:midpoint)
        Iμ = interval(μc, rμ; format=:midpoint)

        θmax_refined = sup(Intθ_g(Iκ, Iμ)[2])
        θmin_refined = inf(Intθ_g(Iκ, Iμ)[1])

        # ------------------------
        # Case 0: no brackets
        # ------------------------
        if isempty(br)
            # certify “no root” by proving g has fixed sign on the whole θ-interval
            FI = eval_subdiv(g, interval(θmin_refined, θmax_refined), Iκ, Iμ; ntiles=50)

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
            root = find_zero(θ -> g(θ, κc, μc), br[1]; atol=1e-12)

            ok_root, data = validate_root_sign_mono(
                g,
                gθ,
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
                    g,
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
            # Case 2: two or more brackets 
            # ------------------------
        else
            root1 = find_zero(θ -> g(θ, κc, μc), br[1]; atol=1e-12)
            root2 = find_zero(θ -> g(θ, κc, μc), br[2]; atol=1e-12)

            # order roots by absolute value
            if abs(root2) < abs(root1)
                (root1, root2) = (root2, root1)
            end

            ok_root1, data1 = validate_root_sign_mono(
                g,
                gθ,
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
                g,
                gθ,
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
                    g,
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


κ_lo, κ_hi = 0.0, 0.5
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 4*2464, 4*89
boxes5_validated, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)

### Second region 
κ_lo, κ_hi = 0.5, 1.5
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 4 * 2464, 4 * 89
boxes5_validated2, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)

### Third region 
κ_lo, κ_hi = 1.5, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 4*2464, 4*89
boxes5_validated3, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


append!(boxes5_validated, boxes5_validated2)
append!(boxes5_validated, boxes5_validated3)

validated_boxes5 = filter(bx -> bx[1] !== :unknown, boxes5_validated)
if isempty(validated_boxes5)
    println("boxes5_validated: no validated boxes found")
else
    bx_min_μ = argmin(bx -> bx[3] - bx[5], validated_boxes5)
    root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx_min_μ
    println("boxes5_validated: smallest validated μ lower edge = $(μc - rμ)")
    println("  label=$(root_lbl), center=(κ=$(κc), μ=$(μc))")
    println("  κ interval=[$(κc - rκ), $(κc + rκ)], μ interval=[$(μc - rμ), $(μc + rμ)]")
end




function plot_regions_all_labels_from_boxes5(boxes;
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    resolution::Int=2000,
    alpha::Real=1.0,
    title::AbstractString="",
    axis_box::Bool=true,

   
    overlay_rect=nothing,                    
    overlay_nx::Int=0,
    overlay_ny::Int=0,
    overlay_color=RGBAf(0.97, 0.98, 1.00, 1),
    overlay_line_color=RGBAf(1, 1, 1, 1),
    overlay_line_alpha::Real=0.95,
    overlay_line_lw::Real=0.9,
    overlay_bg_alpha::Real=0.30,
)
    κlo, κhi = κlims
    μlo, μhi = μlims

    @assert κhi > κlo
    @assert μhi > μlo
    @assert resolution ≥ 2
    @assert (label_type === :root) || (label_type === :stability)
    @assert root_index == 1 || root_index == 2

    κgrid = collect(range(float(κlo), float(κhi); length=resolution))
    μgrid = collect(range(float(μlo), float(μhi); length=resolution))
    dκ = (κhi - κlo) / (resolution - 1)
    dμ = (μhi - μlo) / (resolution - 1)

    @inline clampi(i) = ifelse(i < 1, 1, ifelse(i > resolution, resolution, i))


    @inline function rootcount_color_fixed(lbl::Symbol)
        return (
            lbl === :unknown ? "#1A1A1A" :
            lbl === :zero_roots ? "#F2F2F2" :
            lbl === :one_root ? rootcount_color(:two_roots) :
            lbl === :two_roots ? rootcount_color(:four_roots) :
            rootcount_color(lbl)
        )
    end

    colorfn = (label_type === :root) ? rootcount_color_fixed : stability_color

    
    labels = Symbol[]
    for bx in boxes
        root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        if label_type === :root
            push!(labels, root_lbl)
        else
            if root_lbl === :zero_roots || (root_index == 2 && root_lbl === :one_root)
                # Grey means that the selected root does not exist.
                push!(labels, :zero_roots)
            else
                if root_index == 1
                    ok1 || continue
                    push!(labels, stab1)
                else
                    ok2 || continue
                    push!(labels, stab2)
                end
            end
        end
    end
    labels = sort!(unique(labels); by=String)

    
    masks = Dict{Symbol,Matrix{UInt8}}()
    for lbl in labels
        masks[lbl] = fill(UInt8(0), resolution, resolution)
    end

    
    for bx in boxes
        root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        lbl = if label_type === :root
            root_lbl
        else
            if root_lbl === :zero_roots || (root_index == 2 && root_lbl === :one_root)
                # Reuse the zero-root grey for an absent selected root.
                :zero_roots
            else
                if root_index == 1
                    ok1 || continue
                    stab1
                else
                    ok2 || continue
                    stab2
                end
            end
        end

        κa = max(κc - rκ, κlo)
        κb = min(κc + rκ, κhi)
        μa = max(μc - rμ, μlo)
        μb = min(μc + rμ, μhi)

        i1 = clampi(floor(Int, (κa - κlo) / dκ) + 1)
        i2 = clampi(ceil(Int, (κb - κlo) / dκ) + 1)
        j1 = clampi(floor(Int, (μa - μlo) / dμ) + 1)
        j2 = clampi(ceil(Int, (μb - μlo) / dμ) + 1)

        M = masks[lbl]
        @inbounds for i in i1:i2, j in j1:j2
            M[i, j] = UInt8(1)
        end
    end

    
    fig = Figure()
    ax = Axis(fig[1, 1];
        backgroundcolor=Makie.to_color("#F2F2F2"),
        xlabel=L"\mathbf{\kappa}",
        ylabel=L"\mathbf{\mu}",
        title=title,
        titlesize=AXIS_TITLE_SIZE,
        xlabelsize=AXIS_LABEL_SIZE,
        ylabelsize=AXIS_LABEL_SIZE,
        xticklabelsize=AXIS_TICK_SIZE,
        yticklabelsize=AXIS_TICK_SIZE,
        topspinevisible=true,
        rightspinevisible=true,
        xgridvisible=false,
        ygridvisible=false,
        xtickalign=1,
        ytickalign=1,
        xticksmirrored=true,
        yticksmirrored=true,
        spinewidth=2,
        xtickwidth=2,
        ytickwidth=2,
    )

    xlims!(ax, κlo, κhi)
    ylims!(ax, μlo, μhi)

    labels_plot = [filter(!=(:unknown), labels)...]
    :unknown in labels && push!(labels_plot, :unknown)

    for lbl in labels_plot
        Z = Float32.(masks[lbl])
        col = (lbl === :unknown) ? RGBf(0, 0, 0) : Makie.to_color(colorfn(lbl))
        α = (lbl === :unknown) ? 1.0 : alpha

        heatmap!(ax, κgrid, μgrid, Z;
            colormap=[RGBAf(0, 0, 0, 0), col],
            colorrange=(0f0, 1f0),
            interpolate=false,
            alpha=α,
        )
    end

    
    if overlay_rect !== nothing
        κa, κb, μa, μb = overlay_rect
        κa = max(float(κa), κlo)
        κb = min(float(κb), κhi)
        μa = max(float(μa), μlo)
        μb = min(float(μb), μhi)

        if κb > κa && μb > μa
            poly!(ax,
                Point2f[(κa, μa), (κb, μa), (κb, μb), (κa, μb)];
                color=Makie.to_color(overlay_color),
                alpha=overlay_bg_alpha,
                strokewidth=0,
            )

            if overlay_nx > 0
                Δκg = (κb - κa) / overlay_nx
                for m in 1:(overlay_nx-1)
                    κx = κa + m * Δκg
                    lines!(ax, [κx, κx], [μa, μb];
                        color=overlay_line_color,
                        alpha=overlay_line_alpha,
                        linewidth=overlay_line_lw,
                    )
                end
            end

            if overlay_ny > 0
                Δμg = (μb - μa) / overlay_ny
                for m in 1:(overlay_ny-1)
                    μy = μa + m * Δμg
                    lines!(ax, [κa, κb], [μy, μy];
                        color=overlay_line_color,
                        alpha=overlay_line_alpha,
                        linewidth=overlay_line_lw,
                    )
                end
            end

            lines!(ax, [κa, κa], [μa, μb];
                color=overlay_line_color,
                linewidth=1.4,
                overdraw=true,
            )
            lines!(ax, [κb, κb], [μa, μb];
                color=overlay_line_color,
                linewidth=1.4,
                overdraw=true,
            )
            lines!(ax, [κa, κb], [μa, μa];
                color=overlay_line_color,
                linewidth=1.4,
                overdraw=true,
            )
            lines!(ax, [κa, κb], [μb, μb];
                color=overlay_line_color,
                linewidth=1.4,
                overdraw=true,
            )
        end
    end

    return fig, ax
end


fig1, ax1 = plot_regions_all_labels_from_boxes5(boxes5_validated;
    label_type=:root,
)

save(joinpath(@__DIR__, "..", "figures", "validated", "g_root_regions.png"), fig1)



fig2, ax2 = plot_regions_all_labels_from_boxes5(boxes5_validated;
    label_type=:stability,
    root_index=1
)

save(joinpath(@__DIR__, "..", "figures", "validated", "g_stability_root1.png"), fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes5(boxes5_validated;
    label_type=:stability,
    root_index=2
)

save(joinpath(@__DIR__, "..", "figures", "validated", "g_stability_root2.png"), fig3)


