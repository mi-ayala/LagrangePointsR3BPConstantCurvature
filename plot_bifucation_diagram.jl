include("numerics_f1.jl")
include("numerics_f2.jl")
include("numerics_f3.jl")
include("numerics_f4.jl")


using CairoMakie
function waterfall_regions_from_boxes_clean(boxes_list;
    names=nothing,
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    ms::Real=3.0,
    alpha::Real=0.05,          # translucent
    yspacing::Real=1.0,
    title::AbstractString="Waterfall regions",
)

    κlo, κhi = κlims
    μlo, μhi = μlims
    nsets = length(boxes_list)

    if names === nothing
        names = ["f$i" for i in 1:nsets]
    end

    colorfn = (label_type === :root) ? rootcount_color : stability_color

    fig = Figure()
    ax = Axis3(fig[1, 1],
        title=title,
        xlabel=L"\kappa",
        ylabel="",
        zlabel=L"\mu",
    )

    xlims!(ax, κlo, κhi)
    zlims!(ax, μlo, μhi)
    ylims!(ax, 0, (nsets - 1) * yspacing)

    for (i, boxes) in enumerate(boxes_list)
        y0 = (i - 1) * yspacing

        pts_by_lbl = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()

        for bx in boxes
            root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx

            (κc < κlo || κc > κhi || μc < μlo || μc > μhi) && continue

            if label_type === :root
                # skip zero-root slices completely
                root_lbl === :zero_roots && continue
                lbl = root_lbl
            else
                # stability: skip if this root does not exist
                θ = (root_index == 1) ? θ1 : θ2
                isfinite(θ) || continue
                lbl = (root_index == 1) ? stab1 : stab2
            end

            xy = get!(pts_by_lbl, lbl) do
                (Float64[], Float64[])
            end
            push!(xy[1], κc)
            push!(xy[2], μc)
        end

        for (lbl, (κvals, μvals)) in pts_by_lbl
            isempty(κvals) && continue

            col = Makie.to_color(colorfn(lbl))

            scatter!(ax,
                κvals,
                fill(y0, length(κvals)),
                μvals;
                markersize=ms,
                color=col,
                alpha=alpha,
                transparency=true,
            )
        end
    end

    hidespines!(ax)
    # ax.gridvisible = false

    return fig, ax
end

# figW, axW = waterfall_regions_from_boxes_clean(
#     [boxes1, boxes2, boxes3, boxes4];
#     names=["f1","f2","f3","f4"],
#     label_type=:root,
#     title="",
# )

# save("figures/numerics/waterfall_root_regions.png", figW)

using CairoMakie

# -------------------------
# Existence switches (ALL)
# -------------------------
# boxes tuple format assumed:
# (root_lbl, κc, μc, θ1, θ2, stab1, stab2)
function existence_switch_points(boxes; include_unknown::Bool=false)
    exists(lbl::Symbol) = (lbl !== :zero_roots) && (include_unknown || lbl !== :unknown)

    rows = Dict{Float64,Vector{Tuple{Float64,Bool}}}()
    for bx in boxes
        root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx
        μf = float(μc)
        push!(get!(rows, μf) do
                Tuple{Float64,Bool}[]
            end, (float(κc), exists(root_lbl)))
    end

    # collect all switches as (κ_switch, μ)
    pts = Tuple{Float64,Float64}[]
    for (μ, vec) in rows
        sort!(vec; by=first)  # by κ
        for i in 1:(length(vec)-1)
            κ1, e1 = vec[i]
            κ2, e2 = vec[i+1]
            e1 == e2 && continue
            push!(pts, (0.5 * (κ1 + κ2), μ))
        end
    end

    sort!(pts; by=last)  # by μ
    return pts
end

using CairoMakie

# boxes tuple format assumed:
# (root_lbl::Symbol, κc::Float64, μc::Float64, θ1::Float64, θ2::Float64, stab1::Symbol, stab2::Symbol)

# Return all existence↔no-existence switch points as (κ_switch, μc).
# For each fixed μ, we sort by κ and record every adjacent change of "exists".
function existence_switch_points(boxes; include_unknown::Bool=false)
    exists(lbl::Symbol) = (lbl !== :zero_roots) && (include_unknown || lbl !== :unknown)

    # group κ-values by μ
    byμ = Dict{Float64,Vector{Tuple{Float64,Bool}}}()
    for bx in boxes
        root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx
        μ = float(μc)
        v = get!(byμ, μ) do
            Tuple{Float64,Bool}[]
        end
        push!(v, (float(κc), exists(root_lbl)))
    end

    pts = Tuple{Float64,Float64}[]
    for (μ, vec) in byμ
        sort!(vec; by=first)  # sort by κ
        for i in 1:(length(vec)-1)
            κ1, e1 = vec[i]
            κ2, e2 = vec[i+1]
            e1 == e2 && continue
            push!(pts, (0.5 * (κ1 + κ2), μ))  # midpoint is a simple "switch location"
        end
    end

    sort!(pts; by=last)  # by μ
    return pts
end

# Convert switch points into multiple curves.
# For each μ, there may be 0,1,2,... switches (κ values). Sort them.
# Curve j is the j-th switch across μ.
function switch_points_to_curves(points)
    isempty(points) && return Tuple{Vector{Float64},Vector{Float64}}[]

    byμ = Dict{Float64,Vector{Float64}}()
    for (κ, μ) in points
        v = get!(byμ, μ) do
            Float64[]
        end
        push!(v, κ)
    end

    μs = sort!(collect(keys(byμ)))
    for μ in μs
        sort!(byμ[μ])  # sort κ switches
    end

    mmax = maximum(length(byμ[μ]) for μ in μs)

    curves = Tuple{Vector{Float64},Vector{Float64}}[]
    for j in 1:mmax
        κs = Float64[]
        μs_used = Float64[]
        for μ in μs
            v = byμ[μ]
            length(v) >= j || continue
            push!(κs, v[j])
            push!(μs_used, μ)
        end
        # keep even small curves; caller can decide to skip
        push!(curves, (κs, μs_used))
    end

    return curves
end

using CairoMakie

using CairoMakie

function waterfall_from_curves(curves_list;
    names=nothing,
    yspacing::Real=1.0,
    alpha::Real=0.5,
    line_alpha::Real=0.75,
    linewidth::Real=2.0,
    zbase::Union{Nothing,Real}=nothing,
    colors=nothing,
    title::AbstractString="Waterfall",
    show_labels::Bool=true,
)
    nsets = length(curves_list)
    nsets ≥ 1 || error("curves_list is empty")

    # names handling (avoid BoundsError)
    if names === nothing
        names = ["slice $i" for i in 1:nsets]
    elseif length(names) < nsets
        # pad if user gave too few names
        append!(names, ["slice $i" for i in (length(names)+1):nsets])
    end

    cols = colors === nothing ? Makie.wong_colors() : colors

    # pick baseline in z
    μmins = Float64[]
    for curve in curves_list
        isempty(curve) && continue
        push!(μmins, minimum(last.(curve)))
    end
    z0 = zbase === nothing ? (isempty(μmins) ? 0.0 : minimum(μmins)) : float(zbase)

    fig = Figure()
    ax = Axis3(fig[1, 1],
        title=title,
        xlabel=L"\kappa",
        ylabel="",
        zlabel=L"\mu",
    )

    ax.azimuth[] = 0.5      # front view (along +y)
    ax.elevation[] = 0.0      # no tilt
    # draw from back to front
    for i in nsets:-1:1
        curve = curves_list[i]
        isempty(curve) && continue

        y0 = Float32((i - 1) * yspacing)
        col = cols[mod1(i, length(cols))]

        κ = Float32.(first.(curve))
        μ = Float32.(last.(curve))

        # sort by κ so bands/lines don't zig-zag
        p = sortperm(κ)
        κ = κ[p]
        μ = μ[p]

        y = fill(y0, length(κ))

        points = Point3f.(κ, y, μ)
        base = Point3f.(κ, y, fill(Float32(z0), length(κ)))

        band!(ax, base, points; color=col, alpha=float(alpha), transparency=true)
        lines!(ax, points; color=col, alpha=float(line_alpha), linewidth=float(linewidth))

        if show_labels
            # label near the start of the curve
            text!(ax, names[i];
                position=Point3f(κ[1], y0, μ[end]),
                align=(:left, :bottom),
                fontsize=18
            )
        end
    end

    return fig, ax
end

μconst = 0.0886765
κs = range(0.0, (0.5 * pi)^2; length=2000)

curve = [(float(κ), μconst) for κ in κs]  # Vector{Tuple{Float64,Float64}}
function thin_curve(curve, k)
    return curve[1:k:end]
end

curves_list = [
    thin_curve(curve, 300), thin_curve(curve2, 300), thin_curve(curve3, 300), thin_curve(curve, 300)  # slice 1 has 1 curve
]

fig, ax = waterfall_from_curves(curves_list;
)

save("figures/numerics/waterfall_existence_boundary.png", fig)