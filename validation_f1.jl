using IntervalArithmetic
using Printf
using Roots
using CairoMakie

include("functions_helpers.jl")
include("functions_stability.jl")
include("functions_validation.jl")
include("functions_f1.jl")


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

    boxes = Tuple{Bool,Bool,Float64,Float64,Float64,Float64,Symbol,Symbol}[]
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
        seeded += 1 # Seeded means we found a root 

        # 2) refine the root 
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

        if ok_root

            _, root_label = finalize_root_label(complements_nonzero(f₁, interval(0.0, θmax_refined), data.J, Iκ, Iμ))

            _, stability_label = classify_equilibrium(data.Iexist, Iκ, Iμ)

            if root_label === :unknown
                push!(boxes, (true, true, κc, μc, rκ, rμ, root_label, :unknown))

            else
                push!(boxes, (true, true, κc, μc, rκ, rμ, root_label, stability_label))
                validated += 1
            end


            # println(@sprintf("Validated cell (κ=%.6f, μ=%.6f): root=%s, stability=%s", κc, μc, root_label, stability_label))

        elseif keep_invalid
            push!(boxes, (false, false, κc, μc, rκ, rμ, :unknown, :unknown))
            # println(@sprintf("Failed to validate cell (κ=%.6f, μ=%.6f)", κc, μc))
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

nx, ny = 2*2464, 2*89
# nx, ny = 2500, 100
boxes, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


# ### Plotting

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


using CairoMakie


# --- simple heuristic: choose raster resolution from original interval grid ---
# Oversample > 1 makes smoother boundaries than the interval grid.
@inline function _estimate_raster_res_from_grid(nx::Int, ny::Int;
    oversample::Real=1.5,
    minκ::Int=800, maxκ::Int=4000,
    minμ::Int=250, maxμ::Int=1500,
)
    resκ = clamp(round(Int, oversample * nx), minκ, maxκ)
    resμ = clamp(round(Int, oversample * ny), minμ, maxμ)
    return (resκ, resμ)
end

using CairoMakie

"""
    plot_regions_all_labels_from_boxes(boxes;
        label_type=:root,              # :root or :stability
        κlims=(0.0, (0.5*pi)^2),
        μlims=(0.0, 0.0886765),
        resolution=2000,
        alpha=0.7,
        show_boundary=true,
        boundary_lw=2.0,
        title="",
        axis_box=true)

Overlay (rasterized) regions for *all labels present* in `boxes`, using the
color functions `rootcount_color(::Symbol)` or `stability_color(::Symbol)`.

Boxes are tuples:
    (ok1, ok2, κc, μc, rκ, rμ, root_label::Symbol, stability_label::Symbol)

Special case:
  - label == :unknown is plotted in deep black with alpha = 1.0.
"""
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

    # --- collect labels present ---
    labels = Symbol[]
    for bx in boxes
        ok1, ok2, κc, μc, rκ, rμ, root_lbl, stab_lbl = bx
        push!(labels, (label_type === :root) ? root_lbl : stab_lbl)
    end
    labels = sort!(unique(labels); by=String)

    # --- one mask per label ---
    masks = Dict{Symbol,Matrix{UInt8}}()
    for lbl in labels
        masks[lbl] = fill(UInt8(0), resolution, resolution)
    end

    # --- rasterize each label's union ---
    for bx in boxes
        ok1, ok2, κc, μc, rκ, rμ, root_lbl, stab_lbl = bx
        lbl = (label_type === :root) ? root_lbl : stab_lbl

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

    # --- figure/axis ---
    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel="κ",
        ylabel="μ",
        title=title,
        # xticklabelsize=5,
        # yticklabelsize=5,
        topspinevisible=true,
        rightspinevisible=true,
        xgridvisible=false,
        ygridvisible=false,
        xtickalign=1,
        ytickalign=1,
        xticksmirrored=true,
        yticksmirrored=true,
        xticksize=5,
        yticksize=5,
        spinewidth=1,
        xtickwidth=1,
        ytickwidth=1,
    )

    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end

    xlims!(ax, κlo, κhi)
    ylims!(ax, μlo, μhi)

    # choose color function
    colorfn = (label_type === :root) ? rootcount_color : stability_color

    # IMPORTANT: draw :unknown last so it sits on top, and fully opaque black
    labels_plot = [filter(!=(:unknown), labels)...]
    if :unknown in labels
        push!(labels_plot, :unknown)
    end

    for lbl in labels_plot
        M = masks[lbl]
        Z = Float32.(M)

        # deep black + fully opaque for unknown
        col = (lbl === :unknown) ? RGBf(0, 0, 0) : Makie.to_color(colorfn(lbl))
        α = (lbl === :unknown) ? 1.0 : float(alpha)

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

    return fig, ax
end

fig1, ax = plot_regions_all_labels_from_boxes(boxes; label_type=:root, title="Root-count regions")
# CairoMakie.save("root_regions.pdf", fig)

fig2, ax = plot_regions_all_labels_from_boxes(boxes; label_type=:stability, title="Stability regions")
# CairoMakie.save("stab_regions.pdf", fig)

CairoMakie.save("root_regions.png", fig1)  # or 1
CairoMakie.save("stab_regions.png", fig2)  # or 1


using CairoMakie
using LaTeXStrings

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

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel = L"\mathbf{\kappa}",
        ylabel = L"\mathbf{\mu}",
        title  = title,
        ylabelsize = 20,
        xlabelsize = 20,
        # xticklabelsize = 25,
        # yticklabelsize = 25,
        topspinevisible = true,
        rightspinevisible = true,
        xgridvisible = false,
        ygridvisible = false,
        xtickalign = 1,
        ytickalign = 1,
        xticksmirrored = true,
        yticksmirrored = true,
        # xticksize = 15,
        # yticksize = 15,
        # spinewidth = 4,
        # xtickwidth = 4,
        # ytickwidth = 4,
    )

    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end

    xlims!(ax, κlo, κhi)
    ylims!(ax, μlo, μhi)

    colorfn = (label_type === :root) ? rootcount_color : stability_color

    # draw unknown last so it sits on top
    labels_plot = [filter(!=(:unknown), labels)...]
    if :unknown in labels
        push!(labels_plot, :unknown)
    end

    for lbl in labels_plot
        M = masks[lbl]
        Z = Float32.(M)

        # FORCE unknown to black + fully opaque
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

fig, ax = plot_regions_all_labels_from_boxes(
    boxes;
    label_type = :root,
    κlims = (0.0, (0.5*pi)^2),
    μlims = (0.0, 0.0886765),
    resolution = 1200,        # good balance: clean but light
    alpha = 0.7,
    show_boundary = true,
    boundary_lw = 2.0,
    # title = "Root-count regions",
)

# save("root_regions.pdf", fig)   # vector-friendly + crisp boundaries
# or
save("root_regions.png", fig)