using CairoMakie
using LaTeXStrings


@inline function rootcount_color(lbl::Symbol)
    return (
        lbl === :multiple_roots ? "#08306B" :  # navy
        lbl === :two_roots ? "#08306B" :  # deep blue
        lbl === :one_root ? "#2171B5" :  # medium blue
        lbl === :zero_roots ? "#BDBDBD" :  # medium-light gray (not too light)
        lbl === :unknown ? "#000000" :  # black
        "#000000"
    )
end


# ---------- Stability colors (stable → unstable) ----------
@inline function stability_color(lbl::Symbol)
    return (
        lbl === :zero_roots ? "#BDBDBD" :  # same gray as existence
        lbl === :lyapunov_stable ? "#009E73" :  # teal-green (stable)
        lbl === :elliptic ? "#0072B2" :  # blue (neutral/center-like)
        lbl === :linearly_unstable ? "#E69F00" :  # orange
        lbl === :center_saddle_unstable ? "#E69F00" :  # deeper amber (clear step)
        lbl === :complex_saddle_unstable ? "#D55E00" :  # red-ish (distinct from orange)
        lbl === :saddle_hyperbolic ? "#D55E00" :  # red-orange (most unstable)
        lbl === :unknown ? "#1A1A1A" :  # charcoal
        "#000000"
    )
end

"""
    plot_regions_all_labels_from_boxes(boxes;
        label_type=:root,                 # :root or :stability
        root_index=1,                     # only used when label_type==:stability (1 or 2)
        κlims=(0.0, (0.5*pi)^2),
        μlims=(0.0, 0.0886765),
        resolution=2000,
        alpha=0.7,
        show_boundary=true,
        boundary_lw=2.0,
        title="",
        axis_box=true,
        unknown_box_lw=0.8)

Plots regions from a *regular box cover* by rasterizing the union of boxes per label.

NEW box structure expected:
    (rootcount::Symbol, κc::Float64, μc::Float64, rκ::Float64, rμ::Float64,
     ok1::Bool, ok2::Bool, stab1::Symbol, stab2::Symbol)

Behavior:
- label_type == :root:
    rasterizes regions for each rootcount label.
    (:unknown is NOT rasterized; its boxes are drawn at the end as black rectangles.)
- label_type == :stability:
    for boxes with rootcount == :zero_roots → rasterize as :zero_roots (light gray)
    else use stab1/stab2 depending on root_index.
    (:unknown stability is NOT rasterized; its boxes are drawn at the end as black rectangles.)

Returns: (fig, ax)
Call/signature and return are unchanged.
"""
function plot_regions_all_labels_from_boxes(boxes;
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    resolution::Int=2000,
    alpha::Real=1.0,
    show_boundary::Bool=true,
    boundary_lw::Real=2.0,
    title::AbstractString="",
    axis_box::Bool=true,
    unknown_box_lw::Real=0.8,
)

    κlo, κhi = κlims
    μlo, μhi = μlims
    @assert κhi > κlo
    @assert μhi > μlo
    @assert resolution ≥ 2
    @assert (label_type === :root) || (label_type === :stability)
    @assert (root_index == 1) || (root_index == 2)

    # grid endpoints (Makie likes explicit vectors for heatmap/contour)
    κgrid = collect(range(float(κlo), float(κhi); length=resolution))
    μgrid = collect(range(float(μlo), float(μhi); length=resolution))
    dκ = (κhi - κlo) / (resolution - 1)
    dμ = (μhi - μlo) / (resolution - 1)

    @inline clampi(i) = ifelse(i < 1, 1, ifelse(i > resolution, resolution, i))

    # ------------------- collect labels present (excluding :unknown; we draw unknown as boxes) -------------------
    labels = Symbol[]
    for bx in boxes
        rootcount, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        if label_type === :root
            rootcount === :unknown && continue
            push!(labels, rootcount)
        else
            if rootcount === :zero_roots
                push!(labels, :zero_roots)
            else
                stab = (root_index == 1) ? stab1 : stab2
                stab === :unknown && continue
                push!(labels, stab)
            end
        end
    end
    labels = sort!(unique(labels); by=String)

    # ------------------- allocate one mask per label -------------------
    masks = Dict{Symbol,Matrix{UInt8}}()
    for lbl in labels
        masks[lbl] = fill(UInt8(0), resolution, resolution)
    end

    # ------------------- rasterize each label's union -------------------
    for bx in boxes
        rootcount, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        # choose label for this box depending on mode
        # lbl::Symbol
        if label_type === :root
            rootcount === :unknown && continue
            lbl = rootcount
        else
            if rootcount === :zero_roots
                lbl = :zero_roots
            else
                stab = (root_index == 1) ? stab1 : stab2
                stab === :unknown && continue
                lbl = stab
                # optional gating (comment/uncomment as desired):
                # (root_index == 1 ? ok1 : ok2) || continue
            end
        end

        haskey(masks, lbl) || continue  # safety

        κa = κc - rκ
        κb = κc + rκ
        μa = μc - rμ
        μb = μc + rμ

        (κb < κlo || κa > κhi || μb < μlo || μa > μhi) && continue
        κa = max(κa, κlo)
        κb = min(κb, κhi)
        μa = max(μa, μlo)
        μb = min(μb, μhi)

        i1 = clampi(floor(Int, (κa - κlo) / dκ) + 1)
        i2 = clampi(ceil(Int, (κb - κlo) / dκ) + 1)
        j1 = clampi(floor(Int, (μa - μlo) / dμ) + 1)
        j2 = clampi(ceil(Int, (μb - μlo) / dμ) + 1)

        M = masks[lbl]
        @inbounds for i in i1:i2, j in j1:j2
            M[i, j] = UInt8(1)
        end
    end

    # ------------------- figure/axis -------------------
    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel=L"\mathbf{\kappa}",
        ylabel=L"\mathbf{\mu}",
        title=title,
        ylabelsize=20,
        xlabelsize=20,
        topspinevisible=true,
        rightspinevisible=true,
        xgridvisible=false,
        ygridvisible=false,
        xtickalign=1,
        ytickalign=1,
        xticksmirrored=true,
        yticksmirrored=true,
    )

    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end

    # xlims!(ax, κlo, κhi)
    # ylims!(ax, μlo, μhi)

    # choose color function (you provide these)
    colorfn = (label_type === :root) ? rootcount_color : stability_color

    # ------------------- draw masks (rasterized unions) -------------------
    for lbl in labels
        M = masks[lbl]
        Z = Float32.(M)
        col = Makie.to_color(colorfn(lbl))
        α = float(alpha)

        heatmap!(ax, κgrid, μgrid, Z;
            colormap=[RGBAf(0, 0, 0, 0), col],
            colorrange=(0f0, 1f0),
            interpolate=false,
            alpha=α,
        )

        if show_boundary
            contour!(ax, κgrid, μgrid, Z;
                levels=[0.5f0],
                color=col,
                linewidth=boundary_lw,
                labels=false,
            )
        end
    end

    # ------------------- draw UNKNOWN as actual boxes (black), on top -------------------
    # (This is what you asked: do NOT rasterize unknown pixels; draw the unknown boxes.)
    black = RGBf(0, 0, 0)

    for bx in boxes
        rootcount, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        is_unknown = if label_type === :root
            rootcount === :unknown
        else
            # only unknown stability where there is a root;
            # keep :zero_roots rasterized as gray.
            if rootcount === :zero_roots
                false
            else
                ((root_index == 1) ? stab1 : stab2) === :unknown
            end
        end
        is_unknown || continue

        κa = κc - rκ
        κb = κc + rκ
        μa = μc - rμ
        μb = μc + rμ

        (κb < κlo || κa > κhi || μb < μlo || μa > μhi) && continue
        κa = max(κa, κlo)
        κb = min(κb, κhi)
        μa = max(μa, μlo)
        μb = min(μb, μhi)

        # rectangle outline only (keeps file light compared to filled polys)
        lines!(ax,
            [κa, κb, κb, κa, κa],
            [μa, μa, μb, μb, μa];
            color=black,
            linewidth=unknown_box_lw
        )
    end

    return fig, ax
end

# """
#     plot_regions_all_labels_from_boxes_quick(boxes;
#         label_type=:root,                    # :root or :stability
#         root_index::Int=1,                   # used only when label_type=:stability
#         κlims=(0.0, (0.5*pi)^2),
#         μlims=(0.0, 0.0886765),
#         resolution=2000,
#         alpha=0.7,
#         title="",
#         axis_box=true)

# QUICK version: rasterization ONLY (no boundaries, no per-box overlays).
# Plots all labels present in `boxes` in one CairoMakie axis using heatmap layers.

# Expected box structure:
#     (rootcount::Symbol, κc::Float64, μc::Float64, rκ::Float64, rμ::Float64,
#      ok1::Bool, ok2::Bool, stab1::Symbol, stab2::Symbol)

# - If `label_type == :root`: uses `rootcount` as the label.
# - If `label_type == :stability`: uses `stab{root_index}` as the label,
#   but also rasterizes `rootcount == :zero_roots` as `:zero_roots` (N/A stability),
#   using `stability_color(:zero_roots)`.

# Unknown is forced to deep black and fully opaque (alpha=1).
# Returns: (fig, ax)
# """
function plot_regions_all_labels_from_boxes_quick(boxes;
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    resolution::Int=2000,
    alpha::Real=1.0,
    title::AbstractString="",
    axis_box::Bool=true,
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

    # --- collect labels present ---
    labels = Symbol[]
    for bx in boxes
        root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx
        if label_type === :root
            push!(labels, root_lbl)
        else
            # for stability plots: if no root, label as :zero_roots (N/A stability region)
            if root_lbl === :zero_roots
                push!(labels, :zero_roots)
            else
                push!(labels, root_index == 1 ? stab1 : stab2)
            end
        end
    end
    labels = sort!(unique(labels); by=String)

    # --- one mask per label ---
    masks = Dict{Symbol,Matrix{UInt8}}()
    for lbl in labels
        masks[lbl] = fill(UInt8(0), resolution, resolution)
    end

    # --- rasterize each label's union ---
    for bx in boxes
        root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        lbl = if label_type === :root
            root_lbl
        else
            root_lbl === :zero_roots ? :zero_roots : (root_index == 1 ? stab1 : stab2)
        end

        κa = κc - rκ
        κb = κc + rκ
        μa = μc - rμ
        μb = μc + rμ

        (κb < κlo || κa > κhi || μb < μlo || μa > μhi) && continue
        κa = max(κa, κlo)
        κb = min(κb, κhi)
        μa = max(μa, μlo)
        μb = min(μb, μhi)

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
        xlabel=L"\mathbf{\kappa}",
        ylabel=L"\mathbf{\mu}",
        title=title,
        ylabelsize=20,
        xlabelsize=20,
        topspinevisible=true,
        rightspinevisible=true,
        xgridvisible=false,
        ygridvisible=false,
        xtickalign=1,
        ytickalign=1,
        xticksmirrored=true,
        yticksmirrored=true,
    )
    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end
    xlims!(ax, κlo, κhi)
    ylims!(ax, μlo, μhi)

    colorfn = (label_type === :root) ? rootcount_color : stability_color

    # draw unknown last so it sits on top (still rasterized here)
    labels_plot = [filter(!=(:unknown), labels)...]
    if :unknown in labels
        push!(labels_plot, :unknown)
    end

    for lbl in labels_plot
        M = masks[lbl]
        Z = Float32.(M)

        col = (lbl === :unknown) ? RGBf(0, 0, 0) : Makie.to_color(colorfn(lbl))
        α = (lbl === :unknown) ? 1.0 : float(alpha)

        heatmap!(ax, κgrid, μgrid, Z;
            colormap=[RGBAf(0, 0, 0, 0), col],
            colorrange=(0f0, 1f0),
            interpolate=false,
            alpha=α,
        )
    end

    return fig, ax
end

# function plot_regions_all_labels_from_boxes_quick(boxes;
#     label_type::Symbol=:root,
#     root_index::Int=1,
#     κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
#     μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
#     resolution::Int=2000,
#     alpha::Real=1.0,
#     title::AbstractString="",
#     axis_box::Bool=true,
#     # --- NEW: optional overlay rectangle grid ---
#     overlay_rect=nothing,          # (κa, κb, μa, μb) or nothing
#     overlay_nx::Int=10,            # number of vertical grid lines (incl. edges)
#     overlay_ny::Int=10,            # number of horizontal grid lines (incl. edges)
#     overlay_color=RGBAf(1, 1, 1, 1),  # white by default
#     overlay_lw::Real=1.0,
#     overlay_alpha::Real=0.55,
# )
#     κlo, κhi = κlims
#     μlo, μhi = μlims
#     @assert κhi > κlo
#     @assert μhi > μlo
#     @assert resolution ≥ 2
#     @assert (label_type === :root) || (label_type === :stability)
#     @assert root_index == 1 || root_index == 2

#     κgrid = collect(range(float(κlo), float(κhi); length=resolution))
#     μgrid = collect(range(float(μlo), float(μhi); length=resolution))

#     @inline clampi(i) = ifelse(i < 1, 1, ifelse(i > resolution, resolution, i))

#     # --- collect labels present ---
#     labels = Symbol[]
#     for bx in boxes
#         root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

#         if label_type === :root
#             push!(labels, root_lbl)
#         else
#             # stability plot:
#             # - always allow :zero_roots (N/A stability region)
#             # - otherwise only include stability labels where the chosen root EXISTS
#             if root_lbl === :zero_roots
#                 push!(labels, :zero_roots)
#             else
#                 hasroot = (root_index == 1) ? ok1 : ok2
#                 hasroot || continue
#                 push!(labels, root_index == 1 ? stab1 : stab2)
#             end
#         end
#     end
#     labels = sort!(unique(labels); by=String)

#     # --- one mask per label ---
#     masks = Dict{Symbol,Matrix{UInt8}}()
#     for lbl in labels
#         masks[lbl] = fill(UInt8(0), resolution, resolution)
#     end

#     # --- rasterize each label's union ---
#     for bx in boxes
#         root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

#         # decide label, or skip (transparent) for "root doesn't exist"
#         lbl = if label_type === :root
#             root_lbl
#         else
#             if root_lbl === :zero_roots
#                 :zero_roots
#             else
#                 hasroot = (root_index == 1) ? ok1 : ok2
#                 hasroot || continue          # <-- KEY: don't paint black; just skip
#                 (root_index == 1) ? stab1 : stab2
#             end
#         end

#         κa = κc - rκ
#         κb = κc + rκ
#         μa = μc - rμ
#         μb = μc + rμ

#         (κb < κlo || κa > κhi || μb < μlo || μa > μhi) && continue
#         κa = max(κa, κlo)
#         κb = min(κb, κhi)
#         μa = max(μa, μlo)
#         μb = min(μb, μhi)

#         i1 = clampi(searchsortedfirst(κgrid, κa))
#         i2 = clampi(searchsortedlast( κgrid, κb))
#         j1 = clampi(searchsortedfirst(μgrid, μa))
#         j2 = clampi(searchsortedlast( μgrid, μb))

#         (i1 > i2 || j1 > j2) && continue

#         M = masks[lbl]
#         @inbounds for i in i1:i2, j in j1:j2
#             M[i, j] = UInt8(1)
#         end
#     end

#     fig = Figure()
#     ax = Axis(fig[1, 1];
#         xlabel=L"\mathbf{\kappa}",
#         ylabel=L"\mathbf{\mu}",
#         title=title,
#         ylabelsize=20,
#         xlabelsize=20,
#         topspinevisible=true,
#         rightspinevisible=true,
#         xgridvisible=false,
#         ygridvisible=false,
#         xtickalign=1,
#         ytickalign=1,
#         xticksmirrored=true,
#         yticksmirrored=true,
#     )
#     if axis_box
#         ax.topspinevisible = true
#         ax.rightspinevisible = true
#     end
#     xlims!(ax, κlo, κhi)
#     ylims!(ax, μlo, μhi)

#     colorfn = (label_type === :root) ? rootcount_color : stability_color

#     # draw unknown last so it sits on top
#     labels_plot = [filter(!=(:unknown), labels)...]
#     if :unknown in labels
#         push!(labels_plot, :unknown)
#     end

#     for lbl in labels_plot
#         M = masks[lbl]
#         Z = Float32.(M)

#         col = (lbl === :unknown) ? RGBf(0, 0, 0) : Makie.to_color(colorfn(lbl))
#         α = (lbl === :unknown) ? 1.0 : float(alpha)

#         heatmap!(ax, κgrid, μgrid, Z;
#             colormap=[RGBAf(0, 0, 0, 0), col],
#             colorrange=(0f0, 1f0),
#             interpolate=false,
#             alpha=α,
#         )
#     end

#     # --- overlay a white grid inside a rectangle ---
#     if overlay_rect !== nothing
#         κa, κb, μa, μb = overlay_rect

#         κa = max(float(κa), float(κlo))
#         κb = min(float(κb), float(κhi))
#         μa = max(float(μa), float(μlo))
#         μb = min(float(μb), float(μhi))

#         if κb > κa && μb > μa
#             nx = max(2, overlay_nx)
#             ny = max(2, overlay_ny)

#             col_grid = Makie.to_color(overlay_color)
#             αg = float(overlay_alpha)

#             for κx in range(κa, κb; length=nx)
#                 lines!(ax, [κx, κx], [μa, μb]; color=col_grid, linewidth=overlay_lw, alpha=αg)
#             end
#             for μy in range(μa, μb; length=ny)
#                 lines!(ax, [κa, κb], [μy, μy]; color=col_grid, linewidth=overlay_lw, alpha=αg)
#             end
#         end
#     end

#     return fig, ax
# end
