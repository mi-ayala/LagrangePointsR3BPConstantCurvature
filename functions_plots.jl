using CairoMakie
using LaTeXStrings


# ---------- Root-count colors ----------
@inline function rootcount_color(lbl::Symbol)
    return (
        lbl === :two_roots ? "#2C5AA0" :   # strong blue
        lbl === :one_root ? "#8FBCE6" :   # light blue
        lbl === :zero_roots ? "#D0D0D0" :   # light gray
        lbl === :unknown ? "#000000" :   # black
        "#000000"    # fallback
    )
end


# ---------- Stability colors (stable → unstable) ----------
@inline function stability_color(lbl::Symbol)
    return (
        lbl === :lyapunov_stable ? "#1B9E77" :  # green (stable)
        lbl === :elliptic ? "#A6D854" :  # yellow-green
        lbl === :linearly_unstable ? "#E67E22" :  # orange
        lbl === :saddle_hyperbolic ? "#ED1B2F" :  # red (most unstable)
        lbl === :unknown ? "#000000" :  # black
        "#000000"
    )
end




function plot_regions_all_labels_from_boxes(boxes;
    label_type::Symbol=:root,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    resolution::Int=2000,
    alpha::Real=0.7,
    show_boundary::Bool=true,
    boundary_lw::Real=2.0,
    title::AbstractString="",
    axis_box::Bool=true,
)

    κlo, κhi = κlims
    μlo, μhi = μlims
    @assert κhi > κlo
    @assert μhi > μlo
    @assert resolution ≥ 2
    @assert (label_type === :root) || (label_type === :stability)

    κgrid = collect(range(float(κlo), float(κhi); length=resolution))
    μgrid = collect(range(float(μlo), float(μhi); length=resolution))
    dκ = (κhi - κlo) / (resolution - 1)
    dμ = (μhi - μlo) / (resolution - 1)

    @inline clampi(i) = ifelse(i < 1, 1, ifelse(i > resolution, resolution, i))

    labels = Symbol[]
    for bx in boxes
        ok1, ok2, κc, μc, rκ, rμ, root_lbl, stab_lbl = bx
        push!(labels, (label_type === :root) ? root_lbl : stab_lbl)
    end
    labels = sort!(unique(labels); by=String)

    masks = Dict{Symbol,Matrix{UInt8}}()
    for lbl in labels
        masks[lbl] = fill(UInt8(0), resolution, resolution)
    end

    for bx in boxes
        ok1, ok2, κc, μc, rκ, rμ, root_lbl, stab_lbl = bx
        lbl = (label_type === :root) ? root_lbl : stab_lbl

        κa = κc - rκ; κb = κc + rκ
        μa = μc - rμ; μb = μc + rμ

        (κb < κlo || κa > κhi || μb < μlo || μa > μhi) && continue
        κa = max(κa, κlo); κb = min(κb, κhi)
        μa = max(μa, μlo); μb = min(μb, μhi)

        i1 = clampi(floor(Int, (κa - κlo) / dκ) + 1)
        i2 = clampi(ceil(Int,  (κb - κlo) / dκ) + 1)
        j1 = clampi(floor(Int, (μa - μlo) / dμ) + 1)
        j2 = clampi(ceil(Int,  (μb - μlo) / dμ) + 1)

        M = masks[lbl]
        @inbounds for i in i1:i2, j in j1:j2
            M[i, j] = UInt8(1)
        end
    end

    # One plot sized to fit 2 columns × 3 rows per page
    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel = L"\mathbf{\kappa}",
        ylabel = L"\mathbf{\mu}",
        title  = title,
        ylabelsize = 24,
        xlabelsize = 24,
        xticklabelsize = 18,
        yticklabelsize = 18,
        topspinevisible = true,
        rightspinevisible = true,
        xgridvisible = false,
        ygridvisible = false,
        xtickalign = 1,
        ytickalign = 1,
        xticksmirrored = true,
        yticksmirrored = true,
        xticksize = 10,
        yticksize = 10,
        spinewidth = 2,
        xtickwidth = 2,
        ytickwidth = 2,
    )

    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end

    # xlims!(ax, κlo, κhi)
    # ylims!(ax, μlo, μhi)

    colorfn = (label_type === :root) ? rootcount_color : stability_color

    labels_plot = [filter(!=(:unknown), labels)...]
    if :unknown in labels
        push!(labels_plot, :unknown)
    end

    for lbl in labels_plot
        M = masks[lbl]
        Z = Float32.(M)

        col = (lbl === :unknown) ? RGBf(0,0,0) : Makie.to_color(colorfn(lbl))
        α   = (lbl === :unknown) ? 1.0 : float(alpha)

        heatmap!(ax, κgrid, μgrid, Z;
            colormap = [RGBAf(0,0,0,0), col],
            colorrange = (0f0, 1f0),
            interpolate = false,
            alpha = α,
        )

        if show_boundary
            contour!(ax, κgrid, μgrid, Z;
                levels = [0.5f0],
                color = col,
                linewidth = boundary_lw,
                labels = false,
            )
        end
    end

    return fig, ax
end