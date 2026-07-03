
using CairoMakie
using LaTeXStrings
using RC3BP

const AXIS_TICK_SIZE = RC3BP.AXIS_TICK_SIZE
const AXIS_TITLE_SIZE = RC3BP.AXIS_TITLE_SIZE
const GRID_LINE_WIDTH = RC3BP.GRID_LINE_WIDTH
const GRID_LINE_COLOR = RC3BP.GRID_LINE_COLOR

include("numerics_f5.jl")
include("numerics_f5neg.jl")
append!(boxes5, boxes5neg)

function plot_mu_slice_theta_over_sqrtkappa_triangular(boxes, μ;
    use::Symbol=:stability,
    μtol::Real=Inf,
    κlims::Tuple{<:Real,<:Real}=(-1.0, 1.0),
    alpha::Real=1.0,
    ms::Real=4.0,
    title::AbstractString="",
    axis_box::Bool=true,
    switch_tol::Real=0.1,
    max_switches::Int=2,
    yellow_lbl::Union{Symbol,Nothing}=:center_saddle_unstable,
)
    @assert use === :stability || use === :rootcount

    κlo, κhi = κlims

    # --- pick closest μ slice ---
    μvals = Float64[]
    for bx in boxes
        _, _, μc, _, _, _, _ = bx
        push!(μvals, μc)
    end

    @assert !isempty(μvals) "boxes is empty"

    idx = argmin(abs.(μvals .- float(μ)))
    μstar = μvals[idx]

    if abs(μstar - float(μ)) > float(μtol)
        error("No μ slice within μtol=$(μtol). Closest is μstar=$(μstar).")
    end

    # --- figure/axis ---
    fig = Figure()

    ax = Axis(fig[1, 1];
        xlabel=L"\mathbf{\kappa}",
        titlesize=AXIS_TITLE_SIZE,
        xticklabelsize=AXIS_TICK_SIZE,
        yticklabelsize=AXIS_TICK_SIZE,
        topspinevisible=true,
        rightspinevisible=true,
        xgridvisible=true,
        ygridvisible=true,
        xgridstyle=:dot,
        ygridstyle=:dot,
        xgridwidth=GRID_LINE_WIDTH,
        ygridwidth=GRID_LINE_WIDTH,
        xgridcolor=GRID_LINE_COLOR,
        ygridcolor=GRID_LINE_COLOR,
        xtickalign=1,
        ytickalign=1,
        xticksmirrored=true,
        yticksmirrored=true,
        spinewidth=2,
        xtickwidth=2,
        ytickwidth=2,
    )

    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end

    xlims!(ax, κlo, κhi)

    ybottom, ytop = -0.1, 3.5
    ylims!(ax, ybottom, ytop)

    # group points by label, one scatter per label
    pts_by_lbl = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()

    # all plotted points, independent of root ordering
    # entries are (κ, y, label)
    all_pts = Tuple{Float64,Float64,Symbol}[]

    function _add_point!(lbl::Symbol, κc::Float64, y::Float64)
        xy = get!(pts_by_lbl, lbl) do
            (Float64[], Float64[])
        end

        push!(xy[1], κc)
        push!(xy[2], y)

        push!(all_pts, (κc, y, lbl))

        return nothing
    end

    # (:two_roots, κc, μc, root1, root2, stability_root1, stability_root2)
    for bx in boxes
        root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx

        abs(μc - μstar) > 1e-9 && continue
        (κc < κlo || κc > κhi) && continue
        κc == 0 && continue

        invsqrtκ = inv(sqrt(abs(κc)))

        # root 1
        if isfinite(θ1)
            y1 = θ1 * invsqrtκ
            lbl1 = (use === :stability) ? stab1 : root_lbl
            _add_point!(lbl1, κc, y1)
        end

        # root 2
        if isfinite(θ2)
            y2 = θ2 * invsqrtκ
            lbl2 = (use === :stability) ? stab2 : root_lbl
            _add_point!(lbl2, κc, y2)
        end
    end

    # --- plot points by label ---
    labels = collect(keys(pts_by_lbl))
    labels_plot = [filter(!=(:unknown), labels)...]

    if :unknown in labels
        push!(labels_plot, :unknown)
    end

    println("Available labels: ", labels_plot)

    for lbl in labels_plot
        κs, ys = pts_by_lbl[lbl]

        col = if lbl === :unknown
            RGBf(0, 0, 0)
        elseif use === :stability
            Makie.to_color(stability_color(lbl))
        else
            Makie.to_color(rootcount_color(lbl))
        end

        println("Label ", lbl, " has color ", col)

        α = (lbl === :unknown) ? 1.0 : float(alpha)

        scatter!(ax, κs, ys;
            markersize=ms,
            color=col,
            alpha=α,
        )
    end

    # ------------------------------------------------------------
    # Detect closest color switches geometrically.
    #
    # For each close pair of differently labelled points:
    #   - store the two endpoints, so we can draw a connector
    #     and hide the visible gap;
    #   - use the UPPER endpoint for the vertical line.
    # ------------------------------------------------------------

    switch_candidates = Tuple{
        Float64,
        Float64,Float64,
        Float64,Float64,
        Float64,Float64
    }[]
    # entries are:
    # (distance, κ_switch, y_switch, κ1, y1, κ2, y2)

    κscale = κhi - κlo
    yscale = ytop - ybottom

    for i in 1:(length(all_pts) - 1)
        κi, yi, lbli = all_pts[i]

        for j in (i + 1):length(all_pts)
            κj, yj, lblj = all_pts[j]

            lbli === lblj && continue

            d = sqrt(((κi - κj) / κscale)^2 + ((yi - yj) / yscale)^2)

            if d < switch_tol

                # Use the upper point as the switch point.
                # This makes the vertical line end at the higher branch.
                if yi >= yj
                    κs, ys = κi, yi
                else
                    κs, ys = κj, yj
                end

                push!(switch_candidates, (d, κs, ys, κi, yi, κj, yj))
            end
        end
    end

    # closest label changes first
    sort!(switch_candidates, by = p -> p[1])

    clean_switches = Tuple{
        Float64,Float64,
        Float64,Float64,
        Float64,Float64
    }[]
    # entries are:
    # (κ_switch, y_switch, κ1, y1, κ2, y2)

    for (d, κs, ys, κ1, y1, κ2, y2) in switch_candidates
        already_have_this_switch = false

        for (κold, yold, _, _, _, _) in clean_switches
            d_old = sqrt(((κs - κold) / κscale)^2 + ((ys - yold) / yscale)^2)

            if d_old < switch_tol
                already_have_this_switch = true
                break
            end
        end

        if !already_have_this_switch
            push!(clean_switches, (κs, ys, κ1, y1, κ2, y2))
        end

        if length(clean_switches) >= max_switches
            break
        end
    end

    # Draw connectors to hide gaps, then vertical lines from the box floor.
    for (κs, ys, κ1, y1, κ2, y2) in clean_switches

        # connector between the two differently colored nearby points
        lines!(ax,
            [κ1, κ2],
            [y1, y2];
            color=:grey,
            linewidth=2.5,
            linestyle=:dot,
        )

        # vertical line from floor of box to the chosen switch point
        lines!(ax,
            [κs, κs],
            [ybottom, ys];
            color=:grey,
            linewidth=2.5,
            linestyle=:dot,
        )
    end

    # ------------------------------------------------------------
    # Extra vertical line at the highest point of the chosen yellow label.
    #
    # Important:
    # yellow_lbl must be the actual label symbol, not the color name.
    # Check the printed "Available labels" and "Label ... has color ..."
    # lines to identify which label is yellow.
    # ------------------------------------------------------------

    if yellow_lbl !== nothing
        if haskey(pts_by_lbl, yellow_lbl)
            κ_yellow, y_yellow = pts_by_lbl[yellow_lbl]

            idx_yellow_top = argmax(y_yellow)

            κs = κ_yellow[idx_yellow_top]
            ys = y_yellow[idx_yellow_top]

            println("Yellow label found: ", yellow_lbl)
            println("Highest yellow point: κ = ", κs, ", y = ", ys)

            lines!(ax,
                [κs, κs],
                [ybottom, ys];
                color=:grey,
            linewidth=2.5,
            linestyle=:dot,
            )
        else
            println("yellow_lbl = ", yellow_lbl, " was not found in pts_by_lbl.")
            println("Available labels are: ", collect(keys(pts_by_lbl)))
        end
    else
        println("yellow_lbl is nothing, so no extra yellow vertical line was drawn.")
    end

    return fig, ax, μstar
end

κ_lo, κ_hi = -(0.5 * pi)^2, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, .25


# # Plot
fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    label_type=:root,
    ms=4,
    # κlims = (κ_lo, κ_hi),
    μlims = (μ_lo, μ_hi)
)

save(joinpath(@__DIR__, "..", "figures", "numerics", "f5_root_regions_numerics.png"), fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    label_type=:stability,
    root_index=1,
    κlims = (κ_lo, κ_hi),
    μlims = (μ_lo, μ_hi),
    ms=4
)

save(joinpath(@__DIR__, "..", "figures", "numerics", "f5_stability_root1_numerics.png"), fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    label_type=:stability,
    root_index=2,
    κlims = (κ_lo, κ_hi),
    μlims = (μ_lo, μ_hi),
    ms=4
)

save(joinpath(@__DIR__, "..", "figures", "numerics", "f5_stability_root2_numerics.png"), fig3)



function plot_regions_all_labels_from_boxes_numerics_cross(boxes;
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, .25),
    alpha::Real=0.5,
    ms::Real=3.0,
    title::AbstractString="",
    axis_box::Bool=true,
)
    @assert (label_type === :root) || (label_type === :stability)
    @assert root_index == 1 || root_index == 2

    κlo, κhi = κlims
    μlo, μhi = μlims
    @assert κhi > κlo
    @assert μhi > μlo

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel=L"\mathbf{\kappa}",
        ylabel=L"\mathbf{\mu}",
        title=title,
        titlesize=AXIS_TITLE_SIZE,
        ylabelsize=AXIS_LABEL_SIZE,
        xlabelsize=AXIS_LABEL_SIZE,
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

    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end
    xlims!(ax, κlo, κhi)
    ylims!(ax, μlo, μhi)

    colorfn = (label_type === :root) ? rootcount_color : stability_color

    pts_by_lbl = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()

    for bx in boxes
        root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx

        (κc < κlo || κc > κhi || μc < μlo || μc > μhi) && continue

        if label_type === :root
            lbl = root_lbl
        else
            if root_lbl === :zero_roots
                lbl = :zero_roots
            else
                θ = (root_index == 1) ? θ1 : θ2
                isnan(θ) && continue  # skip if that specific root doesn't exist
                lbl = (root_index == 1) ? stab1 : stab2
            end
        end

        xy = get!(pts_by_lbl, lbl) do
            (Float64[], Float64[])
        end
        push!(xy[1], κc)
        push!(xy[2], μc)
    end

    labels = collect(keys(pts_by_lbl))
    labels_plot = [filter(!=(:unknown), labels)...]
    if :unknown in labels
        push!(labels_plot, :unknown)
    end

    for lbl in labels_plot
        κvals, μvals = pts_by_lbl[lbl]

        col = (lbl === :unknown) ? RGBf(0, 0, 0) : Makie.to_color(colorfn(lbl))
        α = (lbl === :unknown) ? 1.0 : float(alpha)

        scatter!(ax, κvals, μvals;
            markersize=ms,
            color=col,
            alpha=α,
        )
    end

    lines!(ax, [κlo, κhi], [0.04];
        color=:black,
        linewidth=2.2,
        # linestyle=:dash
    )

        lines!(ax, [ -0.01], [μlo, μhi];
        color=:black,
        linewidth=2.2,
        # linestyle=:dash
    )


    return fig, ax
end

fig2, ax2 = plot_regions_all_labels_from_boxes_numerics_cross(boxes5;
    κlims=(-(0.5 * pi)^2, (0.5 * pi)^2),
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f5_stability_root1_negative_numerics_together.png", fig2)


# fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa_triangular(boxes5, 0.01;
#     use=:stability,
#     ms=5,
#     κlims=(-1.0, 1.0)
# )

# save(joinpath(@__DIR__, "..", "figures", "numerics", "triangular_bif_plot_mu01.png"), fig3)

# fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa_triangular(boxes5, 0.06;
#     use=:stability,
#     ms=5,
#     κlims=(-1.0, 1.0)
# )

# save(joinpath(@__DIR__, "..", "figures", "numerics", "triangular_bif_plot_mu06.png"), fig3)


