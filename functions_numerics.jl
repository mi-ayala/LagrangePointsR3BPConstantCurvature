"""
    classify_equilibrium_fp(θ, κ, μ; verbose=false)

Floating-point version of the stability classifier.

Returns:
    (label::Symbol, det_val, lam2_val, Δ_val)
"""
function classify_equilibrium_f(θ, κ, μ; verbose::Bool=false)
    det_val = detλ(θ, κ, μ)
    lam2_val = λ₂(θ, κ, μ)
    Δ_val = Δ(θ, κ, μ)

    label =
        det_val < 0 ? :saddle_hyperbolic :
        lam2_val > 0 ? :lyapunov_stable :
        Δ_val > 0 ? :elliptic :
        Δ_val < 0 ? :linearly_unstable :
        :unknown

    if verbose
        println("θ=$θ κ=$κ μ=$μ → label=$label, det=$det_val, λ₂=$lam2_val, Δ=$Δ_val")
    end

    return label
end

function plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type::Symbol = :root,
    root_index::Int = 1,
    κlims::Tuple{<:Real,<:Real} = (0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real} = (0.0, 0.0886765),
    alpha::Real = 0.7,
    ms::Real = 3.0,
    title::AbstractString = "",
    axis_box::Bool = true,
)
    @assert (label_type === :root) || (label_type === :stability)
    @assert root_index == 1 || root_index == 2

    κlo, κhi = κlims
    μlo, μhi = μlims
    @assert κhi > κlo
    @assert μhi > μlo

    # ---- figure/axis (same style) ----
    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel = L"\mathbf{\kappa}",
        ylabel = L"\mathbf{\mu}",
        title  = title,
        ylabelsize = 20,
        xlabelsize = 20,
        topspinevisible = true,
        rightspinevisible = true,
        xgridvisible = false,
        ygridvisible = false,
        xtickalign = 1,
        ytickalign = 1,
        xticksmirrored = true,
        yticksmirrored = true,
    )
    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end
    xlims!(ax, κlo, κhi)
    ylims!(ax, μlo, μhi)

    # color function
    colorfn = (label_type === :root) ? rootcount_color : stability_color

    # ---- group points by label so we do one scatter! per label (fast) ----
    pts_by_lbl = Dict{Symbol, Tuple{Vector{Float64}, Vector{Float64}}}()

    for bx in boxes
        root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx

        # clip to window early
        (κc < κlo || κc > κhi || μc < μlo || μc > μhi) && continue

    
        if label_type === :root
            lbl = root_lbl
        else
            # stability plot: only plot boxes where the chosen root exists
            θ = (root_index == 1) ? θ1 : θ2
            isnan(θ) && continue
            lbl = (root_index == 1) ? stab1 : stab2
        end

        xy = get!(pts_by_lbl, lbl) do
            (Float64[], Float64[])
        end
        push!(xy[1], κc)
        push!(xy[2], μc)
    end

    # draw unknown last so it sits on top
    labels = collect(keys(pts_by_lbl))
    labels_plot = [filter(!=(:unknown), labels)...]
    if :unknown in labels
        push!(labels_plot, :unknown)
    end

    for lbl in labels_plot
        κvals, μvals = pts_by_lbl[lbl]

        col = (lbl === :unknown) ? RGBf(0, 0, 0) : Makie.to_color(colorfn(lbl))
        α   = (lbl === :unknown) ? 1.0 : float(alpha)

        scatter!(ax, κvals, μvals;
            markersize = ms,
            color = col,
            alpha = α,
        )
    end

    return fig, ax
end