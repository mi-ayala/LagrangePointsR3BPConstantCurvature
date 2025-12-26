using IntervalArithmetic
using Printf
using Roots

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

            push!(boxes, (true, true, κc, μc, rκ, rμ, root_label, stability_label))
            validated += 1
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

nx, ny = 2464, 89
# nx, ny = 2500, 100
boxes, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


# ### Plotting

# ============================================================
# Color maps for CairoMakie / Plots
# ============================================================

# ---------- Existence / root-count colors ----------
@inline function rootcount_color(lbl::Symbol)
    return (
        lbl === :two_roots  ? "#2C5AA0" :   # strong blue
        lbl === :one_root   ? "#8FBCE6" :   # light blue
        lbl === :zero_roots ? "#D0D0D0" :   # light gray
        lbl === :unknown    ? "#ED1B2F" :   # red (failed / unknown)
                              "#000000"    # fallback
    )
end


# ---------- Stability colors (stable → unstable) ----------
@inline function stability_color(lbl::Symbol)
    return (
        lbl === :lyapunov_stable   ? "#1B9E77" :  # green
        lbl === :elliptic          ? "#A6D854" :  # yellow-green
        lbl === :linearly_unstable ? "#E67E22" :  # orange
        lbl === :saddle_hyperbolic ? "#B35806" :  # dark orange / brown
        lbl === :unknown           ? "#ED1B2F" :  # red
                                    "#000000"
    )
end


using CairoMakie

"""
    plot_region_from_boxes(boxes;
        label=:one_root,
        label_type=:root,              # :root or :stability
        κlims=(0.0, (0.5*pi)^2),
        μlims=(0.0, 0.0886765),
        resolution=800,
        color=nothing,                 # if nothing, tries rootcount_color/stability_color; else any Makie color
        alpha=0.7,
        show_boundary=true,
        boundary_lw=2.0,
        title="",
        axis_box=true)

Rasterize (for visualization) the UNION of all boxes matching a given label,
then plot as an image + optional boundary contour.

Boxes are assumed to be tuples:
    (ok1, ok2, κc, μc, rκ, rμ, root_label::Symbol, stability_label::Symbol)

Returns: (fig, ax, M, κgrid, μgrid)
where M is a UInt8 mask of size (resolution,resolution) with 0/1 entries.
"""
function plot_region_from_boxes(boxes;
    label::Symbol = :one_root,
    label_type::Symbol = :root,
    κlims::Tuple{<:Real,<:Real} = (0.0, (0.5*pi)^2),
    μlims::Tuple{<:Real,<:Real} = (0.0, 0.0886765),
    resolution::Int = 2000,
    color = nothing,
    alpha::Real = 0.7,
    show_boundary::Bool = true,
    boundary_lw::Real = 2.0,
    title::AbstractString = "",
    axis_box::Bool = true,
)

    κlo, κhi = κlims
    μlo, μhi = μlims
    @assert κhi > κlo
    @assert μhi > μlo
    @assert resolution ≥ 2

    # grids (for contour and for mapping coordinates -> pixel indices)
    κgrid = range(float(κlo), float(κhi); length=resolution)
    μgrid = range(float(μlo), float(μhi); length=resolution)
    dκ = (κhi - κlo) / (resolution - 1)
    dμ = (μhi - μlo) / (resolution - 1)

    # raster mask: M[i,j] corresponds to (κgrid[i], μgrid[j])
    M = fill(UInt8(0), resolution, resolution)

    # helper: clamp index into 1:resolution
    @inline clampi(i) = ifelse(i < 1, 1, ifelse(i > resolution, resolution, i))

    # choose color (only if you didn't pass one)
    col = color
    if col === nothing
        if label_type === :root
            # expects you have rootcount_color(sym)::Color somewhere; fallback if not
            col = @isdefined(rootcount_color) ? rootcount_color(label) : :dodgerblue
        else
            col = @isdefined(stability_color) ? stability_color(label) : :orange
        end
    end

    # rasterize union of boxes with requested label
    for bx in boxes
        # unpack
        ok1, ok2, κc, μc, rκ, rμ, root_lbl, stab_lbl = bx
        lbl = (label_type === :root) ? root_lbl : stab_lbl
        lbl === label || continue

        # if you want to require ok flags, uncomment one of these:
        # ok1 || continue
        # (ok1 && ok2) || continue

        # box extents
        κa = κc - rκ; κb = κc + rκ
        μa = μc - rμ; μb = μc + rμ

        # skip if no overlap with plotting window
        (κb < κlo || κa > κhi || μb < μlo || μa > μhi) && continue

        # clip to window
        κa = max(κa, κlo); κb = min(κb, κhi)
        μa = max(μa, μlo); μb = min(μb, μhi)

        # convert to index ranges (inclusive)
        i1 = clampi(floor(Int, (κa - κlo)/dκ) + 1)
        i2 = clampi(ceil(Int,  (κb - κlo)/dκ) + 1)
        j1 = clampi(floor(Int, (μa - μlo)/dμ) + 1)
        j2 = clampi(ceil(Int,  (μb - μlo)/dμ) + 1)

        @inbounds for i in i1:i2, j in j1:j2
            M[i,j] = UInt8(1)
        end
    end

    # plot
    fig = Figure()
    ax = Axis(fig[1,1];
        xlabel = "κ",
        ylabel = "μ",
        title = title,
    )
    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end

    # IMPORTANT: image! expects endpoints, not StepRangeLen
    image!(ax, (κlo, κhi), (μlo, μhi), M';
        colormap = [:transparent, col],
        interpolate = false,
        alpha = alpha
    )

    if show_boundary
        contour!(ax, κgrid, μgrid, M';
            levels = [0.5],
            color = col,
            linewidth = boundary_lw
        )
    end

    xlims!(ax, κlo, κhi)
    ylims!(ax, μlo, μhi)

    return fig, ax, M, κgrid, μgrid
end




plt= plot_region_from_boxes(
    boxes;
    label = :one_root,
    label_type = :root,
    title = "Existence: one root"
)