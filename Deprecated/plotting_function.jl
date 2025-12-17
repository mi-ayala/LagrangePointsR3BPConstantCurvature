### Plotting function 
"""
    plot_validated_boxes(boxes;
        color=:green, alpha=0.4, edgecolor=:black,
        markersize=4, title="Validated boxes",
        show_bounding_box=true)

Plot validated boxes of the form (ok::Bool, κ, μ, rκ, rμ),
centered at (κ, μ) with half-width rκ and half-height rμ.

Adds the bounding parameter region:
    κ ∈ [0, (0.5π)^2],  μ ∈ [0, μT≈0.0886765].
"""
function plot_validated_boxes(boxes;
    color=:green, alpha=0.4, edgecolor=:black,
    markersize=4, title="Validated boxes",
    show_bounding_box=true)

    μT = 0.0886765
    κ_lo, κ_hi = 0.0, (0.5 * π)^2
    μ_lo, μ_hi = 0.0, μT

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        title=title, legend=false,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    # draw all validated boxes
    for (ok, κ, μ, rκ, rμ) in boxes
        ok || continue
        xs = [κ - rκ, κ + rκ, κ + rκ, κ - rκ, κ - rκ]
        ys = [μ - rμ, μ - rμ, μ + rμ, μ + rμ, μ - rμ]
        Plots.plot!(plt, xs, ys; color=color, alpha=alpha, linecolor=edgecolor, lw=0.5)
        # Plots.scatter!(plt, [κ], [μ]; color=edgecolor, ms=markersize)
    end

    # bounding box
    if show_bounding_box
        Xb = [κ_lo, κ_hi, κ_hi, κ_lo, κ_lo]
        Yb = [μ_lo, μ_lo, μ_hi, μ_hi, μ_lo]
        Plots.plot!(plt, Xb, Yb; color=:transparent, linecolor=:red, lw=3, label="Bounding region")
    end

    Plots.display(plt)
    return plt
end

"""
    plot_validated_boxes(boxes;
        mode=:validated,                        # :validated | :invalid | :both
        color_valid=:green, color_invalid=:orange,
        alpha=0.4, edgecolor=:black, lw=0.5,
        markersize=4, title="Boxes",
        show_bounding_box=true)

Plot boxes of the form (ok::Bool, κ, μ, rκ, rμ).
`mode=:validated` (default) plots only ok=true,
`mode=:invalid` plots only ok=false,
`mode=:both` plots both (different colors).

Adds the bounding region κ ∈ [0, (0.5π)^2], μ ∈ [0, μT≈0.0886765].
"""
function plot_validated_boxes(boxes;
    mode::Symbol=:validated,
    color_valid=:green, color_invalid=:orange,
    alpha::Real=0.4, edgecolor=:black, lw::Real=3,
    markersize::Real=4, title::AbstractString="Boxes",
    show_bounding_box::Bool=true)

    μT = 0.0886765
    κ_lo, κ_hi = 0.0, (0.5 * π)^2
    μ_lo, μ_hi = 0.0, μT

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        title=title, legend=false,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    # helper to draw one box
    draw_box!(κ, μ, rκ, rμ; col) = begin
        xs = [κ - rκ, κ + rκ, κ + rκ, κ - rκ, κ - rκ]
        ys = [μ - rμ, μ - rμ, μ + rμ, μ + rμ, μ - rμ]
        Plots.plot!(plt, xs, ys; color=col, alpha=alpha, linecolor=edgecolor, lw=lw)
        # Plots.scatter!(plt, [κ], [μ]; color=edgecolor, ms=markersize)
    end

    if mode == :validated || mode == :both
        for (ok, κ, μ, rκ, rμ) in boxes
            ok || continue
            draw_box!(κ, μ, rκ, rμ; col=color_valid)
        end
    end
    if mode == :invalid || mode == :both
        for (ok, κ, μ, rκ, rμ) in boxes
            ok && continue
            draw_box!(κ, μ, rκ, rμ; col=color_invalid)
        end
    end

    if show_bounding_box
        Xb = [κ_lo, κ_hi, κ_hi, κ_lo, κ_lo]
        Yb = [μ_lo, μ_lo, μ_hi, μ_hi, μ_lo]
        Plots.plot!(plt, Xb, Yb; color=:transparent, linecolor=:red, lw=3)
    end

    display(plt)
    return plt
end

"""
    plot_validated_boxes(
        boxes;
        mode=:validated,                    # :validated | :invalid | :both
        color_valid=:green, color_invalid=:red,
        alpha=0.4, edgecolor=:black, lw=0.5,
        markersize=3, title="Boxes",
        show_bounding_box=true,
        # NEW:
        y_exaggerate::Real=1.0,             # multiply rμ for drawing only
        min_frac_h::Real=0.0,               # enforce min visual height = min_frac_h*(μ_hi-μ_lo)
        plot_centers_invalid::Bool=true,    # scatter centers of invalid boxes
        plot_centers_valid::Bool=false,     # (optional) centers of valid boxes
        zoom::Union{Nothing,Tuple{Real,Real,Real,Real}}=nothing # (κlo, κhi, μlo, μhi)
    )

Visualize boxes (ok, κ, μ, rκ, rμ). To make tiny boxes visible, set
`y_exaggerate > 1` or `min_frac_h > 0`. These affect drawing only.
"""

"""
    plot_validated_boxes(
        boxes;
        mode=:validated,                    # :validated | :invalid | :both
        color_valid=:green, color_invalid=:red,
        alpha=0.4, edgecolor=:black, lw=0.5,
        markersize=3, title="Boxes",
        show_bounding_box=true,
        # NEW:
        y_exaggerate::Real=1.0,             # multiply rμ for drawing only
        min_frac_h::Real=0.0,               # enforce min visual height = min_frac_h*(μ_hi-μ_lo)
        plot_centers_invalid::Bool=true,    # scatter centers of invalid boxes
        plot_centers_valid::Bool=false,     # (optional) centers of valid boxes
        zoom::Union{Nothing,Tuple{Real,Real,Real,Real}}=nothing # (κlo, κhi, μlo, μhi)
    )

Visualize boxes (ok, κ, μ, rκ, rμ). To make tiny boxes visible, set
`y_exaggerate > 1` or `min_frac_h > 0`. These affect drawing only.
"""
function plot_validated_boxes(boxes;
    mode::Symbol=:validated,
    color_valid=:green, color_invalid=:red,
    alpha::Real=0.4, edgecolor=:black, lw::Real=0.5,
    markersize::Real=1, title::AbstractString="Boxes",
    show_bounding_box::Bool=true,
    y_exaggerate::Real=1.0,
    min_frac_h::Real=0.0,
    plot_centers_invalid::Bool=true,
    plot_centers_valid::Bool=false,
    zoom::Union{Nothing,Tuple{Real,Real,Real,Real}}=nothing)

    μT = 0.0886765
    κ_lo_full, κ_hi_full = 0.0, (0.5 * π)^2
    μ_lo_full, μ_hi_full = 0.0, μT

    # axes limits
    κ_lo, κ_hi, μ_lo, μ_hi = isnothing(zoom) ? (κ_lo_full, κ_hi_full, μ_lo_full, μ_hi_full) : zoom
    μ_span = μ_hi - μ_lo
    min_h = max(0.0, min_frac_h) * μ_span

    plt = Plots.plot(; xlabel="κ", ylabel="μ", title=title, legend=false,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    # helper: draw one rectangle with visual tweaks on height only
    function draw_box!(κ, μ, rκ, rμ; col)
        rμ_vis = max(rμ * y_exaggerate, min_h / 2)   # visual-only height
        xs = [κ - rκ, κ + rκ, κ + rκ, κ - rκ, κ - rκ]
        ys = [μ - rμ_vis, μ - rμ_vis, μ + rμ_vis, μ + rμ_vis, μ - rμ_vis]
        Plots.plot!(plt, xs, ys; color=col, alpha=alpha, linecolor=edgecolor, lw=lw)
    end

    # optionally collect centers for scatter overlays
    κ_invalid = Float64[]
    μ_invalid = Float64[]
    κ_valid = Float64[]
    μ_valid = Float64[]

    # draw according to mode
    if mode == :validated || mode == :both
        for (ok, κ, μ, rκ, rμ) in boxes
            ok || continue
            draw_box!(κ, μ, rκ, rμ; col=color_valid)
            plot_centers_valid && (push!(κ_valid, κ);
            push!(μ_valid, μ))
        end
    end
    if mode == :invalid || mode == :both
        for (ok, κ, μ, rκ, rμ) in boxes
            ok && continue
            draw_box!(κ, μ, rκ, rμ; col=color_invalid)
            plot_centers_invalid && (push!(κ_invalid, κ);
            push!(μ_invalid, μ))
        end
    end

    # overlay centers (helps when rectangles are too thin)
    if plot_centers_invalid && !isempty(κ_invalid)
        Plots.scatter!(plt, κ_invalid, μ_invalid; color=color_invalid, ms=markersize)
    end
    if plot_centers_valid && !isempty(κ_valid)
        Plots.scatter!(plt, κ_valid, μ_valid; color=color_valid, ms=markersize)
    end

    # outer bounding box (full region), unless zooming makes it confusing
    if show_bounding_box
        Xb = [κ_lo_full, κ_hi_full, κ_hi_full, κ_lo_full, κ_lo_full]
        Yb = [μ_lo_full, μ_lo_full, μ_hi_full, μ_hi_full, μ_lo_full]
        Plots.plot!(plt, Xb, Yb; color=:transparent, linecolor=:red, lw=3)
    end

    display(plt)
    return plt
end
