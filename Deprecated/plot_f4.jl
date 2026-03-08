include("numerics_f4.jl")

function plot_mu_slice_theta_over_sqrtkappa(boxes, μ;
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
        xlabel=L"\mathbf{\kappa}",
        ylabel=L"\mathbf{\theta}/\sqrt{\mathbf{\kappa}}",
        title = title == "" ? "μ ≈ $(μstar)" : title,
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
    ylims!(ax, -2.5,2.5)

    # group points by label (one scatter per label)
    pts_by_lbl = Dict{Symbol, Tuple{Vector{Float64}, Vector{Float64}}}()

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

        μc == μstar || continue
        (κc < κlo || κc > κhi) && continue
        κc <= 0 && continue

        invsqrtκ = inv(sqrt(κc))

        # root 1
        if isfinite(θ1)
            y1 = (θ1-2pi) * invsqrtκ
            lbl1 = (use === :stability) ? stab1 : root_lbl
            _add_point!(lbl1, κc, y1)
        end

        # root 2
        if isfinite(θ2)
            y2 = (θ2-2pi) * invsqrtκ
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

    return fig, ax, μstar
end


fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa(boxes4, 0.01;
    use=:stability,
    ms=4,
)

save("figures/numerics/collinearf4_bif_plot_mu01.png", fig3)

fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa(boxes4, 0.06;
    use=:stability,
    ms=4,
)

save("figures/numerics/collinearf4_bif_plot_mu06.png", fig3)

