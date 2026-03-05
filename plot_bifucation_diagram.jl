include("numerics_f1.jl")
include("numerics_f2.jl")
include("numerics_f3.jl")
include("numerics_f4.jl")


"""
    existence_switch_curves_from_boxes(boxes;
        zero_label=:zero_roots,
        include_unknown=false,
        κ_key=2, μ_key=3, rootlbl_key=1,
        tolμ=0.0)

Given `boxes` with entries like (root_lbl, κc, μc, ...),
extract all κ-locations where existence switches along each μ-slice.

- "exists" means root_lbl != zero_label
- "no root" means root_lbl == zero_label
- If include_unknown=false, any :unknown entries are ignored (do not create switches).
- Returns: Vector{Vector{Tuple{Float64,Float64}}} where each inner vector is a curve:
    [(κ, μ), (κ, μ), ...] sorted by μ.

This handles multiple switches per μ by producing multiple curves.
"""
function existence_switch_curves_from_boxes(boxes;
    zero_label::Symbol=:zero_roots,
    include_unknown::Bool=false,
    κ_key::Int=2,
    μ_key::Int=3,
    rootlbl_key::Int=1,
    tolμ::Real=0.0,
)
    isempty(boxes) && return Vector{Vector{Tuple{Float64,Float64}}}()

    # --- collect (μ -> vector of (κ, exists?)) ---
    byμ = Dict{Float64,Vector{Tuple{Float64,Bool}}}()

    for bx in boxes
        root_lbl = bx[rootlbl_key]::Symbol
        κc = float(bx[κ_key])
        μc = float(bx[μ_key])

        if !include_unknown && root_lbl === :unknown
            continue
        end

        exists = (root_lbl !== zero_label) && (root_lbl !== :unknown)
        push!(get!(byμ, μc) do
                Tuple{Float64,Bool}[]
            end, (κc, exists))
    end

    μs = sort!(collect(keys(byμ)))
    isempty(μs) && return Vector{Vector{Tuple{Float64,Float64}}}()

    # --- for each μ slice, compute all switch κ's (midpoints between adjacent κ grid points) ---
    switches_per_μ = Vector{Vector{Float64}}(undef, length(μs))

    for (iμ, μc) in enumerate(μs)
        v = byμ[μc]
        sort!(v, by=first)   # sort by κ
        ks = [p[1] for p in v]
        ex = [p[2] for p in v]

        sw = Float64[]
        for j in 1:length(ks)-1
            if ex[j] != ex[j+1]
                # midpoint in κ between adjacent grid points
                push!(sw, 0.5 * (ks[j] + ks[j+1]))
            end
        end
        switches_per_μ[iμ] = sw
    end

    # --- stitch switches across μ into curves ---
    curves = Vector{Vector{Tuple{Float64,Float64}}}()

    # active curves: each is last κ plus stored points
    active = Vector{Tuple{Float64,Vector{Tuple{Float64,Float64}}}}()

    for (iμ, μc) in enumerate(μs)
        sw = switches_per_μ[iμ]

        # Optionally merge μ-slices that are "the same" within tolμ (rarely needed)
        # (You can ignore tolμ by leaving it 0.0)
        # We still treat each μ separately; tolμ is mainly for future extension.

        if isempty(sw)
            # no switch at this μ: finalize all active curves
            append!(curves, [c[2] for c in active])
            empty!(active)
            continue
        end

        sort!(sw)

        # Greedy matching of switch points to existing curves by closest κ
        used = falses(length(sw))
        new_active = Vector{Tuple{Float64,Vector{Tuple{Float64,Float64}}}}()

        for (κlast, pts) in active
            # find nearest unused switch
            bestj = 0
            bestd = Inf
            for j in eachindex(sw)
                used[j] && continue
                d = abs(sw[j] - κlast)
                if d < bestd
                    bestd = d
                    bestj = j
                end
            end
            if bestj == 0
                # no continuation: finalize
                push!(curves, pts)
            else
                used[bestj] = true
                κnew = sw[bestj]
                push!(pts, (κnew, μc))
                push!(new_active, (κnew, pts))
            end
        end

        # start new curves for any unused switches
        for j in eachindex(sw)
            used[j] && continue
            κnew = sw[j]
            pts = Tuple{Float64,Float64}[(κnew, μc)]
            push!(new_active, (κnew, pts))
        end

        active = new_active
    end

    # finalize remaining
    append!(curves, [c[2] for c in active])

    # sort points in each curve by μ just in case
    for c in curves
        sort!(c, by=last)
    end

    return curves
end
function existence_switch_curve_flat_from_boxes(boxes; kwargs...)
    curves = existence_switch_curves_from_boxes(boxes; kwargs...)
    out = Tuple{Float64,Float64}[]
    for c in curves
        append!(out, c)
    end
    return out
end


"""
    add_top_horizontals(curve; κlo=0.0, κhi, μtol=1e-12, n=200)

Given `curve::Vector{Tuple{Float64,Float64}}` of (κ, μ),
find μmax, then find κL=min κ and κR=max κ among points with μ≈μmax.

Return a new curve where we add two horizontal segments at μ=μmax:
  [κlo → κL] and [κR → κhi].

- `μtol` controls how we detect "points at μmax" (floating tolerance).
- `n` is the number of points used to sample each horizontal segment.
"""
function add_top_horizontals(curve::Vector{Tuple{Float64,Float64}};
    κlo::Float64=0.0,
    κhi::Float64,
    μtol::Float64=1e-12,
    n::Int=2000,
)
    isempty(curve) && error("curve is empty")

    # ensure Float64 and finite
    c = [(Float64(κ), Float64(μ)) for (κ, μ) in curve if isfinite(κ) && isfinite(μ)]
    isempty(c) && error("curve has no finite points")

    μmax = maximum(last.(c))

    # points "at the top"
    top = [(κ, μ) for (κ, μ) in c if abs(μ - μmax) ≤ μtol]
    isempty(top) && error("No points detected at μmax. Increase μtol (try 1e-8).")

    κL = minimum(first.(top))
    κR = maximum(first.(top))

    κlo ≤ κL || error("Need κlo ≤ κL, got κlo=$κlo, κL=$κL")
    κR ≤ κhi || error("Need κR ≤ κhi, got κR=$κR, κhi=$κhi")

    # Build the two horizontal segments (sampled)
    left_seg = [(κ, μmax) for κ in range(κlo, κL; length=max(2, n))]
    right_seg = [(κ, μmax) for κ in range(κR, κhi; length=max(2, n))]

    # Return: left horizontal + original curve + right horizontal
    return vcat(left_seg, c, right_seg)
end

"""
If the curve has a single top maximum μmax, extend it horizontally
from κlo to the first κ where μ = μmax.

curve :: Vector{Tuple{Float64,Float64}}  (κ, μ)
"""
function complete_left_to_single_max(curve; κlo=0.0, n=2000, μtol=1e-12)
    μmax = maximum(last.(curve))

    top = [p for p in curve if abs(p[2] - μmax) ≤ μtol]
    isempty(top) && error("No points detected at μmax (maybe increase μtol).")

    κL = minimum(first.(top))  # leftmost κ on the top

    left = [(κ, μmax) for κ in range(κlo, κL; length=n)]

    return vcat(left, curve)
end


curve2 = existence_switch_curve_flat_from_boxes(boxes2)
curve3 = existence_switch_curve_flat_from_boxes(boxes3)

"""
Resample a curve (κ, μ) to `n` points using piecewise-linear interpolation in κ.

- curve: Vector{Tuple{Float64,Float64}} or any Real pairs
- n: number of output points
- κlo/κhi: optional override for κ range (defaults to min/max κ in curve)
- remove_duplicates: if multiple points share the same κ, keep the last one
Returns: Vector{Tuple{Float64,Float64}} of length n
"""
function resample_curve_linear(curve;
    n::Int=200,
    κlo=nothing,
    κhi=nothing,
    remove_duplicates::Bool=true,
)
    n ≥ 2 || error("n must be ≥ 2")
    isempty(curve) && error("curve is empty")

    # sort by κ
    pts = [(float(p[1]), float(p[2])) for p in curve]
    sort!(pts; by=first)

    # optionally remove duplicate κ (keep last)
    if remove_duplicates
        compact = Tuple{Float64,Float64}[]
        i = 1
        while i ≤ length(pts)
            κ = pts[i][1]
            j = i
            while j ≤ length(pts) && pts[j][1] == κ
                j += 1
            end
            push!(compact, pts[j-1])   # keep last at this κ
            i = j
        end
        pts = compact
    end

    length(pts) ≥ 2 || error("Need at least 2 distinct κ points")

    κmin = pts[1][1]
    κmax = pts[end][1]
    a = κlo === nothing ? κmin : float(κlo)
    b = κhi === nothing ? κmax : float(κhi)
    b > a || error("Need κhi > κlo")

    κgrid = collect(range(a, b; length=n))

    out = Vector{Tuple{Float64,Float64}}(undef, n)

    # pointer walk (O(n + m))
    j = 1
    for i in 1:n
        κ = κgrid[i]

        # clamp outside original range
        if κ ≤ pts[1][1]
            out[i] = (κ, pts[1][2])
            continue
        elseif κ ≥ pts[end][1]
            out[i] = (κ, pts[end][2])
            continue
        end

        while pts[j+1][1] < κ
            j += 1
        end

        κ1, μ1 = pts[j]
        κ2, μ2 = pts[j+1]

        # linear interpolation
        t = (κ - κ1) / (κ2 - κ1)
        μ = (1 - t) * μ1 + t * μ2
        out[i] = (κ, μ)
    end

    return out
end


curve2 = complete_left_to_single_max(curve2)
curve3 = add_top_horizontals(curve3; κlo=0.0, κhi=(0.5 * pi)^2)

curve2 = resample_curve_linear(curve2; n=600)
curve3 = resample_curve_linear(curve3; n=1000)


function waterfall_from_curves(curves_list;
    names=nothing,
    yspacing::Real=1.0,
    alpha::Real=0.2,
    line_alpha::Real=0.95,
    linewidth::Real=2.0,
    zbase::Union{Nothing,Real}=nothing,
    colors=nothing,
    title::AbstractString="",
    show_labels::Bool=false,
)
    nsets = length(curves_list)
    nsets ≥ 1 || error("curves_list is empty")

    # names handling (avoid BoundsError)
    if names === nothing
        names = ["slice $i" for i in 1:nsets]
    elseif length(names) < nsets
        # pad if user gave too few names
        # append!(names, ["slice $i" for i in (length(names)+1):nsets])
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

   # no tilt

    # draw from back to front
    for i in nsets:-1:1
        curve = curves_list[i]
        isempty(curve) && continue

        y0 = Float32((i - 1) * yspacing) + 1
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
κs = range(0.0, (0.5 * pi)^2; length=1000)

curve = [(float(κ), μconst) for κ in κs]  # Vector{Tuple{Float64,Float64}}
function thin_curve(curve, k)
    return curve[1:k:end]
end

curves_list = [
    curve, curve2, curve3, curve  # slice 1 has 1 curve
]

fig, ax = waterfall_from_curves(
    curves_list;
    colors=["#6B7FD7", "#08306B", "#08306B", "#6B7FD7"]
)

save("figures/numerics/waterfall_existence_boundary.png", fig)




using CairoMakie

function waterfall_from_curves2(curves_list;
    names=nothing,
    yspacing::Real=1.0,
    alpha::Real=0.2,
    line_alpha::Real=0.95,
    linewidth::Real=2.0,
    zbase::Union{Nothing,Real}=nothing,
    colors=nothing,
    title::AbstractString="",
    show_labels::Bool=false,
    clean_axes::Bool=true,
    orthographic::Bool=true,
)
    nsets = length(curves_list)
    nsets ≥ 1 || error("curves_list is empty")

    # names
    if names === nothing
        names = ["slice $i" for i in 1:nsets]
    elseif length(names) < nsets
        append!(names, ["slice $i" for i in (length(names)+1):nsets])
    end

    cols = colors === nothing ? Makie.wong_colors() : colors

    # baseline in μ (like minimum(spectrum))
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
        ylabel="",           # y is slice index
        zlabel=L"\mu",
    )

    # “standing at κ label” (front view)
    ax.azimuth = -π/2
    ax.elevation = 0

    # orthographic projection (no perspective distortion)
    # if orthographic
    #     ax.projectiontype = Makie.Orthographic  # NOTE: no parentheses
    # end

    if clean_axes
        # remove grids
        ax.xgridvisible = false
        ax.ygridvisible = false
        ax.zgridvisible = false

        # hide 3D panels (the translucent walls)
        ax.xypanelvisible = false
        ax.yzpanelvisible = false
        ax.xzpanelvisible = false

        # keep spines or hide them (choose: here we keep x/z, hide y)
        ax.xspinesvisible = false
        ax.yspinesvisible = false
        ax.zspinesvisible = false

        # hide the y axis ticks/labels entirely
        ax.yticksvisible = false
        ax.yticklabelsvisible = false
        ax.ylabelvisible = false
        ax.yticks = (Float64[], String[])
    end

    # draw back-to-front
    for i in nsets:-1:1
        curve = curves_list[i]
        isempty(curve) && continue

        y0  = Float32((i - 1) * yspacing) + 1f0
        col = cols[mod1(i, length(cols))]

        κ = Float32.(first.(curve))
        μ = Float32.(last.(curve))

        # sort by κ
        p = sortperm(κ)
        κ = κ[p]
        μ = μ[p]

        y = fill(y0, length(κ))
        points = Point3f.(κ, y, μ)
        base   = Point3f.(κ, y, fill(Float32(z0), length(κ)))

        band!(ax, base, points; color=col, alpha=float(alpha), transparency=true)
        lines!(ax, points; color=col, alpha=float(line_alpha), linewidth=float(linewidth))

        if show_labels
            text!(ax, names[i];
                position=Point3f(κ[1], y0, μ[end]),
                align=(:left, :bottom),
                fontsize=18
            )
        end
    end

    return fig, ax
end

fig, ax = waterfall_from_curves2(
    curves_list;
    colors=["#6B7FD7", "#08306B", "#08306B", "#6B7FD7"]
)

save("figures/numerics/waterfall_existence_boundary2.png", fig)

