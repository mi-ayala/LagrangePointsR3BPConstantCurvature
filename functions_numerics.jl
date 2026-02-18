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
        det_val < 0 ? :center_saddle_unstable :
        lam2_val > 0 ? :lyapunov_stable :
        Δ_val > 0 ? :elliptic :
        Δ_val < 0 ? :complex_saddle_unstable :
        :unknown

    if verbose
        println("θ=$θ κ=$κ μ=$μ → label=$label, det=$det_val, λ₂=$lam2_val, Δ=$Δ_val")
    end

    return label
end


"""
    classify_center_equilibrium_num(x0, κ, μ; verbose=false)

Numeric (Float64) version of the center-equilibrium test:

Λ  = ( sin(√κ q₁) / ( μ sin(√κ q₂) ) )^(1/3)
d₁ = asin( Λ * sin(x0) )
d₂ = x0

Steps:
1) If det(d₁,d₂) < 0  → :center_saddle_unstable
2) else if tr(d₁,d₂) > 0 → :lyapunov_stable
3) else if C(d₁,d₂) > 0 → :elliptic
   else → :complex_saddle_unstable

Returns:
    (ok::Bool, label::Symbol, det::Float64, tr::Float64, C::Float64, d1::Float64, d2::Float64, Λ::Float64)

`ok` is true unless we hit a domain/NaN issue (e.g. asin arg outside [-1,1]).
"""
function classify_center_equilibrium_num(x0::Real, κ::Real, μ::Real; verbose::Bool=false)
    x0 = float(x0)
    κ = float(κ)
    μ = float(μ)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)

    Λ = cbrt(sin(r * q1) / (μ * sin(r * q2)))

    arg = Λ * sin(x0)
    if !(isfinite(arg) && -1.0 ≤ arg ≤ 1.0)
        verbose && println("asin argument not in [-1,1] or not finite: ", arg)
        z = NaN
        return false, :unknown, z, z, z, z, x0, Λ
    end

    d1 = asin(arg)
    d2 = x0

    sr = sin(r)
    sr4 = sr^4

    s_rq1 = sin(r * q1)
    s_rq2 = sin(r * q2)

    pref = (3.0 * s_rq1 * s_rq2) / sr4

    B1 = s_rq2 * cos(d1) + s_rq1 * cos(d2)

    B2 = (1.0 + 2.0 * cos(2.0 * d2)) * s_rq1 * cos(d1) +
         (1.0 + 2.0 * cos(2.0 * d1)) * s_rq2 * cos(d2)

    det = pref * B1 * B2

    tr = (1.0 / (sr^2)) * (
        (1.0 - 4.0 * cos(d2)^2) * s_rq1^2
        -
        6.0 * cos(d1) * cos(d2) * s_rq1 * s_rq2
        +
        (1.0 - 4.0 * cos(d1)^2) * s_rq2^2
    )

    m = -(sinsh = (sin(d1)^2 * sin(d2)^2)) /
        ((cos(r) - cos(d1 + d2)) * (cos(r) - cos(d1 - d2)))

    C = m - 4.0 * det

    verbose && println("det=$det  tr=$tr  C=$C")

    if det < 0
        return :center_saddle_unstable
    elseif tr > 0
        return :lyapunov_stable
    elseif C > 0
        return :elliptic
    else
        return :complex_saddle_unstable
    end
end

function plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    alpha::Real=1.0,
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

    return fig, ax
end

using CairoMakie

"""
    plot_mu_slice_theta_over_sqrtkappa(boxes, μ;
        root_index=1,
        use=:stability,              # :stability or :rootcount
        μtol=Inf,
        κlims=(0.0, (0.5*pi)^2),
        alpha=1.0,
        ms=4.0,
        title="",
        axis_box=true)

Pick μ* = closest μc present in `boxes` to the requested `μ`.
Then plot κ (x-axis) vs θ_root_index / √κ (y-axis) for that root only.

- Skips entries where that root does not exist (θ is NaN / not finite).
- Coloring:
    use=:stability  -> stability_color(stab_root_index)
    use=:rootcount  -> rootcount_color(root_lbl)

Returns (fig, ax, μstar).
"""
function plot_mu_slice_theta_over_sqrtkappa(boxes, μ;
    root_index::Int=1,
    use::Symbol=:stability,
    μtol::Real=Inf,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    alpha::Real=1.0,
    ms::Real=4.0,
    title::AbstractString="",
    axis_box::Bool=true,
)
    @assert root_index == 1 || root_index == 2
    @assert use === :stability || use === :rootcount

    κlo, κhi = κlims
    @assert κhi > κlo

    # --- choose closest μ slice present ---
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

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel=L"\mathbf{\kappa}",
        ylabel=L"\mathbf{\theta}/\sqrt{\mathbf{\kappa}}",
        title = title == "" ? "μ ≈ $(μstar) (root $(root_index))" : title,
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

    # group points by label for fast plotting (1 scatter per label)
    pts_by_lbl = Dict{Symbol, Tuple{Vector{Float64}, Vector{Float64}}}()

    for bx in boxes
        root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx

        μc == μstar || continue
        (κc < κlo || κc > κhi) && continue
        κc <= 0 && continue  # avoid θ/sqrt(κ) blow-up at κ=0

        θ = (root_index == 1) ? θ1 : θ2
        isfinite(θ) || continue   # <-- key: skip where this root doesn't exist

        y = θ / sqrt(κc)

        lbl = if use === :stability
            (root_index == 1) ? stab1 : stab2
        else
            root_lbl
        end

        xy = get!(pts_by_lbl, lbl) do
            (Float64[], Float64[])
        end
        push!(xy[1], κc)
        push!(xy[2], y)
    end

    # plot unknown last
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

using CairoMakie

"""
    plot_theta_over_sqrtkappa_3d_from_boxes(boxes;
        root_index=1,
        use=:stability,                 # :stability or :rootcount
        κlims=(0.0, (0.5*pi)^2),
        μlims=(0.0, 0.0886765),
        ms=3.0,
        alpha=1.0,
        title="",
        axis_box=true,
    )

3D plot of (κ, μ, θ/√κ) for *all* boxes, for the selected root.

Box format expected:
    (root_lbl::Symbol, κc::Float64, μc::Float64, θ1::Float64, θ2::Float64, stab1::Symbol, stab2::Symbol)

Rules:
- If root_lbl == :zero_roots → skip (no θ)
- If θ for chosen root is NaN → skip (that specific root doesn’t exist)
- Colors:
    use=:stability  → stability_color(stab_root)
    use=:rootcount  → rootcount_color(root_lbl)
- Unknown forced to black and drawn last.
Returns (fig, ax).
"""
function plot_theta_over_sqrtkappa_3d_from_boxes(boxes;
    root_index::Int=1,
    use::Symbol=:stability,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5*pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    ms::Real=3.0,
    alpha::Real=1.0,
    title::AbstractString="",
    axis_box::Bool=true,
)
    @assert root_index == 1 || root_index == 2
    @assert use === :stability || use === :rootcount

    κlo, κhi = κlims
    μlo, μhi = μlims

    fig = Figure()
    ax = Axis3(fig[1,1];
        xlabel = L"\mathbf{\kappa}",
        ylabel = L"\mathbf{\mu}",
        zlabel = L"\mathbf{\theta}/\sqrt{\mathbf{\kappa}}",
        title  = title,
    )
    xlims!(ax, float(κlo), float(κhi))
    ylims!(ax, float(μlo), float(μhi))

    # group points by label -> one scatter! per label (fast)
    pts_by_lbl = Dict{Symbol, Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}()

    for bx in boxes
        root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx

        (κc < κlo || κc > κhi || μc < μlo || μc > μhi) && continue
        root_lbl === :zero_roots && continue

        θ = (root_index == 1) ? θ1 : θ2
        isnan(θ) && continue

        z = θ / sqrt(κc)

        lbl = if use === :rootcount
            root_lbl
        else
            (root_index == 1) ? stab1 : stab2
        end

        xyz = get!(pts_by_lbl, lbl) do
            (Float64[], Float64[], Float64[])
        end
        push!(xyz[1], κc)
        push!(xyz[2], μc)
        push!(xyz[3], z)
    end

    colorfn = (use === :rootcount) ? rootcount_color : stability_color

    labels = collect(keys(pts_by_lbl))
    labels_plot = [filter(!=(:unknown), labels)...]
    if :unknown in labels
        push!(labels_plot, :unknown)
    end

    for lbl in labels_plot
        κvals, μvals, zvals = pts_by_lbl[lbl]
        col = (lbl === :unknown) ? RGBf(0,0,0) : Makie.to_color(colorfn(lbl))
        α   = (lbl === :unknown) ? 1.0 : float(alpha)

        scatter!(ax, κvals, μvals, zvals;
            markersize=float(ms),
            color=col,
            alpha=α,
        )
    end

    return fig, ax
end