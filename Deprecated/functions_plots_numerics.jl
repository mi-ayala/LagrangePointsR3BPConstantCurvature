using Printf
using Plots


#### Plotting stability functions 
"""
    stability_color(sym::Symbol)

Map stability labels to colors for plotting.
"""
function stability_color(sym::Symbol)
    if sym === :saddle_hyperbolic
        return :red
    elseif sym === :lyapunov_stable
        return :green
    elseif sym === :elliptic
        return :blue
    elseif sym === :linearly_unstable
        return :orange
    elseif sym === :unknown
        return :gray
    elseif sym === :no_root
        return :black
    else
        return :black
    end
end

"""
    plot_stability_centers(results;
        ms=5,
        title="Stability in (κ, μ) plane",
        show_bounding_box=true)

`results` is a Vector of (κ, μ, θstar, stab::Symbol).

Makes a 2D scatter plot in the (κ, μ) plane with colors by stability.
"""
function plot_stability_centers(results;
    ms::Real=5,
    title::AbstractString="Stability in (κ, μ) plane",
    show_bounding_box::Bool=true)

    κ_vals = getindex.(results, 1)
    μ_vals = getindex.(results, 2)
    θ_vals = getindex.(results, 3)
    stabs = getindex.(results, 4)

    # parameter bounds (adjust if needed)
    μT = 0.0886765
    κ_lo, κ_hi = 0.0, (0.5 * π)^2
    μ_lo, μ_hi = 0.0, μT

    # Group by stability
    groups = Dict{Symbol,Vector{Tuple{Float64,Float64}}}()
    for (κ, μ, θ, stab) in results
        # skip points with NaN root if you like
        if isnan(θ)
            continue
        end
        push!(get!(groups, stab, Tuple{Float64,Float64}[]), (κ, μ))
    end

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        title=title, legend=true,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    for (stab, pts) in groups
        xs = first.(pts)
        ys = last.(pts)
        col = stability_color(stab)
        Plots.scatter!(plt, xs, ys;
            ms=ms, color=col, label=String(stab))
    end

    if show_bounding_box
        Xb = [κ_lo, κ_hi, κ_hi, κ_lo, κ_lo]
        Yb = [μ_lo, μ_lo, μ_hi, μ_hi, μ_lo]
        Plots.plot!(plt, Xb, Yb; color=:transparent, linecolor=:black, lw=2, label=false)
    end

    Plots.display(plt)
    return plt
end



invalid_boxes(boxes) = filter(b -> !b[1], boxes)
function show_invalid_summary(boxes; preview=10)
    bad = invalid_boxes(boxes)
    println("Invalid boxes: ", length(bad))
    for (k, b) in enumerate(bad[1:min(preview, end)])
        ok, κ, μ, rκ, rμ = b
        @printf "%3d) κ=[%.6g, %.6g], μ=[%.6g, %.6g]\n" k (κ - rκ) (κ + rκ) (μ - rμ) (μ + rμ)
    end
    return bad
end

function stability_types(boxes)
    labels = getindex.(boxes, 6)     # grab the 6th element of each tuple
    return unique(labels)           # distinct symbols
end

"""
    plot_stability_centers(boxes;
        only_ok=true,
        show_invalid=false,
        μT=0.0886765,
        zoom=nothing,                # = (κ_lo, κ_hi, μ_lo, μ_hi) or nothing
        markersize=4,
        title="Stability (centers)",
        colors=Dict(
            :saddle_hyperbolic => :red,
            :lyapunov_stable   => :green,
            :elliptic          => :blue,
            :linearly_unstable => :orange,
            :unknown           => :gray
        ),
        labels=nothing)

Plot the centers (κ, μ) of boxes, colored by stability label (6th entry).

`boxes` entries are `(ok::Bool, κ::Float64, μ::Float64, rκ::Float64, rμ::Float64, label::Symbol)`.

- `labels = nothing`  → show all labels
- `labels = :elliptic` → only that label
- `labels = [:elliptic, :lyapunov_stable]` → only those labels
"""
function plot_stability_centers(boxes;
    only_ok=true,
    show_invalid=false,
    μT=0.0886765,
    zoom=nothing,
    markersize=4,
    title="Stability (centers)",
    colors=Dict(
        :saddle_hyperbolic => :red,
        :lyapunov_stable => :green,
        :elliptic => :blue,
        :linearly_unstable => :orange,
        :unknown => :gray
    ),
    labels=nothing)

    # normalize labels filter to a Set{Symbol} or nothing
    labelset = if labels === nothing
        nothing
    elseif labels isa Symbol
        Set{Symbol}((labels,))
    elseif labels isa AbstractVector
        Set(Symbol.(labels))
    else
        error("labels must be nothing, a Symbol, or a vector of Symbols.")
    end

    # full default domain
    κ_lo_full, κ_hi_full = 0.0, (0.5 * π)^2
    μ_lo_full, μ_hi_full = 0.0, μT

    if zoom === nothing
        κ_lo, κ_hi, μ_lo, μ_hi = κ_lo_full, κ_hi_full, μ_lo_full, μ_hi_full
    else
        κ_lo, κ_hi, μ_lo, μ_hi = zoom
    end

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        title=title, legend=true,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    # group centers by label
    bylabel = Dict{Symbol,Vector{Tuple{Float64,Float64}}}()

    for (ok, κ, μ, rκ, rμ, label) in boxes
        if only_ok && !ok
            continue
        end
        if !ok && !show_invalid
            continue
        end
        if labelset !== nothing && !(label in labelset)
            continue
        end
        push!(get!(bylabel, label, Tuple{Float64,Float64}[]), (κ, μ))
    end

    # scatter each label with its own color
    for (lbl, pts) in sort(collect(bylabel); by=x -> String(x[1]))
        κs = [p[1] for p in pts]
        μs = [p[2] for p in pts]
        col = get(colors, lbl, :black)
        Plots.scatter!(plt, κs, μs;
            color=col, ms=markersize, alpha=0.9,
            label=String(lbl))
    end

    Plots.display(plt)
    return plt
end

function plot_invalid_centers(boxes; ms=3, color=:red, title="Invalid boxes (centers)")
    bad = invalid_boxes(boxes)
    κs = getindex.(bad, 2)   # centers
    μs = getindex.(bad, 3)
    μT = 0.0886765
    κ_lo, κ_hi = 0.0, (0.5 * π)^2
    μ_lo, μ_hi = 0.0, μT
    plt = Plots.scatter(κs, μs; ms=ms, color=color, xlabel="κ", ylabel="μ",
        title=title, legend=false, xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    Xb = [κ_lo, κ_hi, κ_hi, κ_lo, κ_lo]
    Yb = [μ_lo, μ_lo, μ_hi, μ_hi, μ_lo]
    Plots.plot!(plt, Xb, Yb; color=:transparent, linecolor=:red, lw=3, label="Bounding region")
end

############################
# Overlay version:  plot_invalid_centers!
############################
function plot_invalid_centers!(plt, boxes;
    ms=3,
    color=:red,
    μT=0.0886765,
    zoom=nothing,
    title=nothing,
    show_bounding_box=true)

    bad = invalid_boxes(boxes)
    isempty(bad) && return plt

    κs = getindex.(bad, 2)   # centers
    μs = getindex.(bad, 3)

    # optionally set limits / title
    κ_lo_full, κ_hi_full = 0.0, (0.5 * π)^2
    μ_lo_full, μ_hi_full = 0.0, μT

    if zoom === nothing
        κ_lo, κ_hi, μ_lo, μ_hi = κ_lo_full, κ_hi_full, μ_lo_full, μ_hi_full
    else
        κ_lo, κ_hi, μ_lo, μ_hi = zoom
        Plots.plot!(plt; xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))
    end

    if title !== nothing
        Plots.plot!(plt; title=title)
    end

    # scatter invalid centers
    Plots.scatter!(plt, κs, μs;
        ms=ms, color=color, xlabel="κ", ylabel="μ",
        legend=false)

    # bounding box if requested
    if show_bounding_box
        Xb = [κ_lo_full, κ_hi_full, κ_hi_full, κ_lo_full, κ_lo_full]
        Yb = [μ_lo_full, μ_lo_full, μ_hi_full, μ_hi_full, μ_lo_full]
        Plots.plot!(plt, Xb, Yb;
            color=:transparent, linecolor=:red, lw=3, label=false)
    end

    return plt
end

############################
# Fresh-plot wrapper: plot_invalid_centers
############################
function plot_invalid_centers(boxes;
    ms=3,
    color=:red,
    μT=0.0886765,
    zoom=nothing,
    title="Invalid boxes (centers)",
    show_bounding_box=true)

    κ_lo_full, κ_hi_full = 0.0, (0.5 * π)^2
    μ_lo_full, μ_hi_full = 0.0, μT

    if zoom === nothing
        κ_lo, κ_hi, μ_lo, μ_hi = κ_lo_full, κ_hi_full, μ_lo_full, μ_hi_full
    else
        κ_lo, κ_hi, μ_lo, μ_hi = zoom
    end

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi),
        title=title, legend=false)

    plot_invalid_centers!(plt, boxes;
        ms=ms, color=color, μT=μT,
        zoom=zoom, title=title,
        show_bounding_box=show_bounding_box)

    Plots.display(plt)
    return plt
end



"""
    plot_stability_centers!(plt, boxes;
        only_ok=true,
        show_invalid=false,
        μT=0.0886765,
        zoom=nothing,                # = (κ_lo, κ_hi, μ_lo, μ_hi) or nothing
        markersize=4,
        title=nothing,               # if given, set title once
        colors=Dict(
            :saddle_hyperbolic => :red,
            :lyapunov_stable   => :green,
            :elliptic          => :blue,
            :linearly_unstable => :orange,
            :unknown           => :gray
        ),
        labels=nothing)

Add the centers (κ, μ) of `boxes` to an existing plot `plt`,
colored by the stability label (6th entry).

`boxes` entries are `(ok::Bool, κ::Float64, μ::Float64, rκ::Float64, rμ::Float64, label::Symbol)`.

- `labels = nothing`  → show all labels
- `labels = :elliptic` → only that label
- `labels = [:elliptic, :lyapunov_stable]` → only those labels
"""
function plot_stability_centers!(plt, boxes;
    only_ok=true,
    show_invalid=false,
    μT=0.0886765,
    zoom=nothing,
    markersize=4,
    title=nothing,
    colors=Dict(
        :saddle_hyperbolic => :red,
        :lyapunov_stable => :green,
        :elliptic => :blue,
        :linearly_unstable => :orange,
        :unknown => :gray
    ),
    labels=nothing)

    # normalize labels filter to a Set{Symbol} or nothing
    labelset = if labels === nothing
        nothing
    elseif labels isa Symbol
        Set{Symbol}((labels,))
    elseif labels isa AbstractVector
        Set(Symbol.(labels))
    else
        error("labels must be nothing, a Symbol, or a vector of Symbols.")
    end

    # optionally set axis limits if zoom is provided
    if zoom !== nothing
        κ_lo, κ_hi, μ_lo, μ_hi = zoom
        Plots.plot!(plt; xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))
    end

    # optionally set title (but don’t overwrite if user wants to manage it)
    if title !== nothing
        Plots.plot!(plt; title=title)
    end

    # group centers by label
    bylabel = Dict{Symbol,Vector{Tuple{Float64,Float64}}}()

    for (ok, κ, μ, rκ, rμ, label) in boxes
        if only_ok && !ok
            continue
        end
        if !ok && !show_invalid
            continue
        end
        if labelset !== nothing && !(label in labelset)
            continue
        end
        push!(get!(bylabel, label, Tuple{Float64,Float64}[]), (κ, μ))
    end

    # scatter per label
    for (lbl, pts) in sort(collect(bylabel); by=x -> String(x[1]))
        κs = [p[1] for p in pts]
        μs = [p[2] for p in pts]
        col = get(colors, lbl, :black)
        Plots.scatter!(plt, κs, μs;
            color=col, ms=markersize, alpha=0.9,
            label=String(lbl))
    end

    return plt
end

"""
    plot_stability_centers(boxes; kwargs...)

Create a new plot and call `plot_stability_centers!` on it.
Returns the plot object.
"""
function plot_stability_centers(boxes; kwargs...)
    # default full domain
    plt = Plots.plot(; xlabel="κ", ylabel="μ", legend=true)
    plot_stability_centers!(plt, boxes; kwargs...)
    Plots.display(plt)
    return plt
end

