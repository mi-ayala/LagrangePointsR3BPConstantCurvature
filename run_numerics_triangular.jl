

include("numerics_f5.jl")
include("numerics_f5neg.jl")
append!(boxes5, boxes5neg)


# include("functions_plots.jl")
# # Plot

fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    label_type=:root,
    ms=4,
)

save("figures/numerics/f5_root_regions_numerics.png", fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f5_stability_root1_numerics.png", fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    label_type=:stability,
    root_index=2,
    ms=4
)

save("figures/numerics/f5_stability_root2_numerics.png", fig3)




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
        # ylabelsize=20,
        # xlabelsize=20,
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

    # ### For f3
    # ax.xticks = (
    #     [0, (pi / 6)^2, (pi / 4)^2, (pi / 3)^2, (pi / 2)^2],
    #     [L"0", L"(\pi/6)^2", L"(\pi/4)^2", L"(\pi/3)^2", L"(\pi/2)^2"]
    # )

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
            # IMPORTANT: keep the "no root" region in stability plots
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

    lines!(ax, [κlo, κhi], [0.040064205622887726, 0.040064205622887726];
        color=:black,
        linewidth=3,
        # linestyle=:dot
    )

    lines!(ax, [0.0, 0.0], [μlo, μhi];
        color=:black,
        linewidth=3,
        # linestyle=:dot
    )
    return fig, ax
end


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    # κlims=(-(0.5 * pi)^2, (0.5 * pi)^2),
    κlims=(-1.0, 1.0),
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f5_stability_root1_negative_numerics_together.png", fig2)


function plot_mu_slice_theta_over_sqrtkappa3(boxes, μ;
    use::Symbol=:stability,
    μtol::Real=Inf,
    κlims::Tuple{<:Real,<:Real}=(-2.6, 2.6),   # <-- allow negatives
    alpha::Real=1.0,
    ms::Real=4.0,
    title::AbstractString="",
    axis_box::Bool=true,
)
    @assert use === :stability || use === :rootcount

    κlo, κhi = κlims
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
        spinewidth=3,
        xtickwidth=2,
        ytickwidth=2,
    )

    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end

    xlims!(ax, κlo, κhi)
    ylims!(ax, -3.5, 3.5)

    # --- group points ---
    pts_by_lbl = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()

    function _add_point!(lbl::Symbol, κc::Float64, y::Float64)
        xy = get!(pts_by_lbl, lbl) do
            (Float64[], Float64[])
        end
        push!(xy[1], κc)
        push!(xy[2], y)
        return nothing
    end

    # --- loop boxes ---
    for bx in boxes
        root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx

        abs(μc - μstar) > 1e-9 && continue
        (κc < κlo || κc > κhi) && continue
        κc == 0 && continue   # <-- only skip zero

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

        scatter!(ax, κs, ys;
            markersize=ms,
            color=col,
            alpha=α,
        )
    end

    return fig, ax, μstar
end

fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa3(boxes5, 0.01;
    use=:stability,
    ms=5,
    κlims=(-1.0, 1.0)
)

save("figures/numerics/triangular_bif_plot_mu01.png", fig3)

fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa3(boxes5, 0.06;
    use=:stability,
    ms=5,
    κlims=(-1.0, 1.0)
)

save("figures/numerics/triangular_bif_plot_mu06.png", fig3)

