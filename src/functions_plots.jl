using CairoMakie
set_theme!(Theme(
    font="Computer Modern",
    fontsize=22
))

const AXIS_LABEL_SIZE = 18
const AXIS_TICK_SIZE = 16
const AXIS_TITLE_SIZE = 18
const GRID_LINE_WIDTH = 1.6
const GRID_LINE_COLOR = RGBAf(0.35, 0.35, 0.35, 0.65)

@inline function rootcount_color(lbl::Symbol)

    pal = Makie.to_colormap(:seaborn_colorblind6)

    return (
        lbl === :unknown ? "#1A1A1A" :
        lbl === :zero_roots ? "#F2F2F2" :
        lbl === :one_root ? pal[1] :
        lbl === :two_roots ? pal[3] :
        lbl === :three_roots ? pal[2] :
        lbl === :four_roots ? pal[6] :
        lbl === :five_roots ? pal[4] :
        lbl === :six_roots ? pal[5] : "#000000"
    )
end

@inline function stability_color(lbl::Symbol)

    pal = Makie.to_colormap(:rainbow)

    return (
        lbl === :unknown ? "#1A1A1A" :
        lbl === :zero_roots ? "#F2F2F2" :
        lbl === :lyapunov_stable ? pal[2] :
        lbl === :elliptic ? pal[3] :
        lbl === :center_saddle_unstable ? pal[4] :
        lbl === :complex_saddle_unstable ? pal[5] : "#000000"
    )
end

function plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
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


    return fig, ax
end

function plot_regions_all_labels_from_total_boxes(boxes;
    κlims::Tuple{<:Real,<:Real}=(-(0.5 * pi)^2, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    alpha::Real=0.1,
    ms::Real=3.0,
    title::AbstractString="",
    axis_box::Bool=true,
)
    κlo, κhi = κlims
    μlo, μhi = μlims
    @assert κhi > κlo
    @assert μhi > μlo

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel=L"\mathbf{\kappa}",
        ylabel=L"\mathbf{\mu}",
        title=title,
        topspinevisible=true,
        rightspinevisible=true,
        xgridvisible=false,
        ygridvisible=false,
        titlesize=AXIS_TITLE_SIZE,
        ylabelsize=AXIS_LABEL_SIZE,
        xlabelsize=AXIS_LABEL_SIZE,
        xticklabelsize=AXIS_TICK_SIZE,
        yticklabelsize=AXIS_TICK_SIZE,
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

    pts_by_lbl = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()

    for bx in boxes
        lbl, κc, μc = bx

        (κc < κlo || κc > κhi || μc < μlo || μc > μhi) && continue

        xy = get!(pts_by_lbl, lbl) do
            (Float64[], Float64[])
        end
        push!(xy[1], κc)
        push!(xy[2], μc)
    end

    # --- add the left region [-((0.5π)^2), 0] as :three_roots ---
    left_col = Makie.to_color(rootcount_color(:three_roots))
    if κlo < 0
        poly!(ax,
            Point2f[
                (κlo, μlo),
                (0.0, μlo),
                (0.0, μhi),
                (κlo, μhi),
            ];
            color=left_col,
            alpha=0.5,
            strokewidth=0,
        )
    end

    labels = collect(keys(pts_by_lbl))
    labels_plot = [filter(!=(:unknown), labels)...]
    if :unknown in labels
        push!(labels_plot, :unknown)
    end

    for lbl in labels_plot
        κvals, μvals = pts_by_lbl[lbl]

        col = (lbl === :unknown) ? RGBf(0, 0, 0) : Makie.to_color(rootcount_color(lbl))
        α = (lbl === :unknown) ? 1.0 : float(alpha)

        scatter!(ax, κvals, μvals;
            markersize=ms,
            color=col,
            alpha=α,
        )
    end

    return fig, ax
end

function combine_rootcount_boxes(boxes_list; unknown_wins::Bool=true)
    isempty(boxes_list) && return Tuple{Symbol,Float64,Float64}[]

    # map labels -> counts
    rootcount(lbl::Symbol) = (
        lbl === :zero_roots ? 0 :
        lbl === :one_root ? 1 :
        lbl === :two_roots ? 2 :
        lbl === :three_roots ? 3 :
        lbl === :four_roots ? 4 :
        lbl === :five_roots ? 5 :
        lbl === :six_roots ? 6 :
        lbl === :unknown ? -1 :
        error("Unknown root-count label: $lbl")
    )

    # map counts -> labels
    countlabel(n::Int) = (
        n == 0 ? :zero_roots :
        n == 1 ? :one_root :
        n == 2 ? :two_roots :
        n == 3 ? :three_roots :
        n == 4 ? :four_roots :
        n == 5 ? :five_roots :
        n == 6 ? :six_roots :
        error("Total root count $n not supported")
    )

    # build a dictionary for each boxes array: (κ,μ) -> label
    dicts = Vector{Dict{Tuple{Float64,Float64},Symbol}}()
    for boxes in boxes_list
        d = Dict{Tuple{Float64,Float64},Symbol}()
        for bx in boxes
            lbl = bx[1]
            κc = float(bx[2])
            μc = float(bx[3])
            d[(κc, μc)] = lbl
        end
        push!(dicts, d)
    end

    # use the first grid as reference
    ref = boxes_list[1]
    out = Tuple{Symbol,Float64,Float64}[]

    for bx in ref
        κc = float(bx[2])
        μc = float(bx[3])
        key = (κc, μc)

        total = 0
        bad = false

        for d in dicts
            haskey(d, key) || error("Grid mismatch: missing point $(key) in one boxes array")
            lbl = d[key]
            n = rootcount(lbl)

            if n == -1
                if unknown_wins
                    bad = true
                    break
                else
                    continue
                end
            end

            total += n
        end

        outlbl = bad ? :unknown : countlabel(total)
        push!(out, (outlbl, κc, μc))
    end

    return out
end



function plot_mu_slice_theta_over_sqrtkappa_collinear(boxes, μ;
    use::Symbol=:stability,
    μtol::Real=Inf,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    alpha::Real=1.0,
    ms::Real=4.0,
    title::AbstractString="",
    axis_box::Bool=true,
    switch_tol::Real=0.01,
    max_switches::Int=4,
)
    @assert use === :stability || use === :rootcount

    κlo, κhi = κlims
    μlo, μhi = 0.0, 0.0886765
    @assert κhi > κlo

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
        ylabelsize=AXIS_LABEL_SIZE,
        xlabelsize=AXIS_LABEL_SIZE,
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

    ybottom, ytop = -3.5, 3.5
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
        κc <= 0 && continue

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

    for lbl in labels_plot
        κs, ys = pts_by_lbl[lbl]

        col = if lbl === :unknown
            RGBf(0, 0, 0)
        elseif use === :stability
            Makie.to_color(stability_color(lbl))
        else
            Makie.to_color(rootcount_color(lbl))
        end

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
    # We compare all pairs of plotted points.
    # If two points have different labels and are very close in the
    # plotted (κ, y) coordinates, this is a candidate label switch.
    #
    # Then we sort all candidates by distance and keep only the closest
    # non-duplicate switches.
    # ------------------------------------------------------------

    switch_candidates = Tuple{Float64,Float64,Float64}[]
    # entries are (distance, κ_switch, y_switch)

    κscale = κhi - κlo
    yscale = ytop - ybottom

    for i in 1:(length(all_pts) - 1)
        κi, yi, lbli = all_pts[i]

        for j in (i + 1):length(all_pts)
            κj, yj, lblj = all_pts[j]

            lbli === lblj && continue

            d = sqrt(((κi - κj) / κscale)^2 + ((yi - yj) / yscale)^2)

            if d < switch_tol
                κs = 0.5 * (κi + κj)
                ys = 0.5 * (yi + yj)

                push!(switch_candidates, (d, κs, ys))
            end
        end
    end

    # Closest label changes first
    sort!(switch_candidates, by = p -> p[1])

    clean_switch_pts = Tuple{Float64,Float64}[]

    for (d, κs, ys) in switch_candidates
        already_have_this_switch = false

        for (κold, yold) in clean_switch_pts
            d_old = sqrt(((κs - κold) / κscale)^2 + ((ys - yold) / yscale)^2)

            if d_old < switch_tol
                already_have_this_switch = true
                break
            end
        end

        if !already_have_this_switch
            push!(clean_switch_pts, (κs, ys))
        end

        if length(clean_switch_pts) >= max_switches
            break
        end
    end

    # Draw lines from the bottom of the box to each switch point.
    for (κs, ys) in clean_switch_pts
        lines!(ax,
            [κs, κs],
            [ybottom, ys];
            color=:grey,
            linestyle=:dot,
            linewidth=2.5,
        )
    end

    # --- overlay q1/sqrt(κ) and q2/sqrt(κ) ---
    κplot = collect(range(max(κlo, 1e-8), κhi; length=800))

    q1vals = [q₁(κ, μstar) for κ in κplot]
    q2vals = [q₂(κ, μstar) for κ in κplot]

    lines!(ax, κplot, q1vals;
        color=:black,
        linewidth=2,
        linestyle=:dash,
    )

    lines!(ax, κplot, q2vals;
        color=:black,
        linewidth=2,
        linestyle=:dash,
    )

    return fig, ax, μstar
end



function plot_regions_all_labels_from_boxes(boxes;
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    resolution::Int=2000,
    alpha::Real=1.0,
    title::AbstractString="",
    axis_box::Bool=true,

    # --- pretty analytical region overlay ---
    overlay_rect=nothing,                    # (κa, κb, μa, μb)
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

    # --- collect labels ---
    labels = Symbol[]
    for bx in boxes
        root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx
        if label_type === :root
            push!(labels, root_lbl)
        else
            if root_lbl === :zero_roots
                push!(labels, :zero_roots)
            else
                push!(labels, root_index == 1 ? stab1 : stab2)
            end
        end
    end
    labels = sort!(unique(labels); by=String)

    # --- masks ---
    masks = Dict{Symbol,Matrix{UInt8}}()
    for lbl in labels
        masks[lbl] = fill(UInt8(0), resolution, resolution)
    end

    # --- rasterize boxes ---
    for bx in boxes
        root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        lbl = if label_type === :root
            root_lbl
        else
            root_lbl === :zero_roots ? :zero_roots :
            (root_index == 1 ? stab1 : stab2)
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

    # --- figure ---
    fig = Figure()
    ax = Axis(fig[1, 1];
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

    colorfn = (label_type === :root) ? rootcount_color : stability_color

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

    # --- analytical region overlay ---
    if overlay_rect !== nothing
        κa, κb, μa, μb = overlay_rect
        κa = max(float(κa), κlo)
        κb = min(float(κb), κhi)
        μa = max(float(μa), μlo)
        μb = min(float(μb), μhi)

        if κb > κa && μb > μa
            # background
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


            # rectangle boundary drawn over the axis
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

function plot_unknown_boxes(boxes;
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    resolution::Int=2000,
    alpha::Real=1.0,
    title::AbstractString="",
    axis_box::Bool=true,

    # --- pretty analytical region overlay ---
    overlay_rect=nothing,                    # (κa, κb, μa, μb)
    overlay_nx::Int=0,
    overlay_ny::Int=0,
    overlay_color=RGBAf(0.97, 0.98, 1.00, 1),
    overlay_line_color=RGBAf(1, 1, 1, 1),
    overlay_line_alpha::Real=1.0,
    overlay_line_lw::Real=1.0,
    overlay_bg_alpha::Real=1.0,
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

    # --- unknown mask only ---
    M = fill(UInt8(0), resolution, resolution)

    # --- rasterize only unknown boxes ---
    for bx in boxes
        root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        lbl = if label_type === :root
            root_lbl
        else
            root_lbl === :zero_roots ? :zero_roots :
            (root_index == 1 ? stab1 : stab2)
        end

        lbl === :unknown || continue

        κa = max(κc - rκ, κlo)
        κb = min(κc + rκ, κhi)
        μa = max(μc - rμ, μlo)
        μb = min(μc + rμ, μhi)

        i1 = clampi(floor(Int, (κa - κlo) / dκ) + 1)
        i2 = clampi(ceil(Int, (κb - κlo) / dκ) + 1)
        j1 = clampi(floor(Int, (μa - μlo) / dμ) + 1)
        j2 = clampi(ceil(Int, (μb - μlo) / dμ) + 1)

        @inbounds for i in i1:i2, j in j1:j2
            M[i, j] = UInt8(1)
        end
    end

    # --- figure ---
    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel=L"\mathbf{\kappa}",
        ylabel=L"\mathbf{\mu}",
        title=title,
        topspinevisible=true,
        rightspinevisible=true,
        ylabelsize=AXIS_LABEL_SIZE,
        xlabelsize=AXIS_LABEL_SIZE,
        xticklabelsize=AXIS_TICK_SIZE,
        yticklabelsize=AXIS_TICK_SIZE,
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

    Z = Float32.(M)
    heatmap!(ax, κgrid, μgrid, Z;
        colormap=[RGBAf(0, 0, 0, 0), RGBf(0, 0, 0)],
        colorrange=(0f0, 1f0),
        interpolate=false,
        alpha=1.0,
    )

    # --- analytical region overlay ---
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
