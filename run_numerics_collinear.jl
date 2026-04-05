using RC3BP

include("numerics_f1.jl")
include("numerics_f2.jl")
include("numerics_f3.jl")
include("numerics_f4.jl")

boxes = Any[]
append!(boxes, boxes1)
append!(boxes, boxes2)
append!(boxes, boxes3)
append!(boxes, boxes4)



### Plots f1
fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes1;
    label_type=:root,
    ms=4
)
save("figures/numerics/f1_root_regions_numerics.png", fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes1;
    label_type=:stability,
    root_index=1,
    ms=4
)
save("figures/numerics/f1_stability_numerics.png", fig2)


### Plots f2

fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes2;
    label_type=:root,
    ms=4,
)

save("figures/numerics/f2_root_regions_numerics.png", fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes2;
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f2_stability_root1_numerics.png", fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes2;
    label_type=:stability,
    root_index=2,
    ms=4
)

save("figures/numerics/f2_stability_root2_numerics.png", fig3)


### Plots f3

fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes3;
    label_type=:root,
    ms=4,
)

save("figures/numerics/f3_root_regions_numerics.png", fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes3;
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f3_stability_root1_numerics.png", fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes3;
    label_type=:stability,
    root_index=2,
    ms=4
)

save("figures/numerics/f3_stability_root2_numerics.png", fig3)


### Plots f4

fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes4;
    label_type=:root,
    ms=4
)

save("figures/numerics/f4_root_regions_numerics.png", fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes4;
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f4_stability_numerics.png", fig2)


#### Function to count roots through all the functions.

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

boxes_total = combine_rootcount_boxes([boxes1, boxes2, boxes3, boxes4])



function plot_regions_all_labels_from_total_boxes2(boxes;
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


fig, ax = plot_regions_all_labels_from_total_boxes2(boxes_total;
    ms=4,
    # title="Total number of roots",
    κlims=(-1.0, (0.5 * pi)^2),
)

save("figures/numerics/total_root_regions.png", fig)



#### Bifurcation diagram
function plot_mu_slice_theta_over_sqrtkappa2(boxes, μ;
    use::Symbol=:stability,
    μtol::Real=Inf,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    alpha::Real=1.0,
    ms::Real=4.0,
    title::AbstractString="",
    axis_box::Bool=true,
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
        # xlabel=L"\mathbf{\kappa}",
        # ylabel=L"\mathbf{\theta}/\sqrt{\mathbf{|\kappa|}}",
        # title=title == "" ? "μ ≈ $(μstar)" : title,
        topspinevisible=true,
        rightspinevisible=true,
        xgridvisible=true,
        ygridvisible=true,
        xgridstyle=:dot,
        ygridstyle=:dot,
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
    ylims!(ax, -3.5, 3.5)


    # group points by label (one scatter per label)
    pts_by_lbl = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()

    # helper: add one root point
    function _add_point!(lbl::Symbol, κc::Float64, y::Float64)
        xy = get!(pts_by_lbl, lbl) do
            (Float64[], Float64[])
        end
        push!(xy[1], κc)
        push!(xy[2], y)
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

        scatter!(ax, κs, ys; markersize=ms, color=col, alpha=α)
    end

    # # --- overlay q1/sqrt(κ) and q2/sqrt(κ) ---
    κplot = collect(range(max(κlo, 1e-8), κhi; length=800))
    μplot = collect(range(μlo, μhi; length=800))

    q1vals = [q₁(κ, μstar) for κ in κplot]
    q2vals = [q₂(κ, μstar) for κ in κplot]

    muval = [(pi / 2)^2 for _ in κplot]

    lines!(ax, κplot, q1vals; color=:black, linewidth=2, linestyle=:dash)
    lines!(ax, κplot, q2vals; color=:black, linewidth=2, linestyle=:dash)

    # vlines!(ax,[(pi/2)^2] ; color=:black, linewidth=2)



    return fig, ax, μstar
end



fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa2(boxes, 0.01;
    use=:stability,
    ms=5,
)

save("figures/numerics/collinear_bif_plot_mu01.png", fig3)



fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa2(boxes, 0.06;
    ms=5,
)

save("figures/numerics/collinear_bif_plot_mu06.png", fig3)
