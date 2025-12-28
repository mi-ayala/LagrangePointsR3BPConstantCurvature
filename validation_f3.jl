using IntervalArithmetic
using Printf
using Roots
using CairoMakie

include("functions_helpers.jl")
include("functions_stability.jl")
include("functions_validation.jl")
include("functions_f3.jl")


function classify_f3_rootcount(Iθm::Interval, Iκ::Interval, Iμ::Interval; verbose::Bool=false)

    Ileft = interval(-pi / 2)
    Iright = -(sqrt(Iκ) * q₁(Iκ, Iμ))
    Iallowed = interval(inf(Ileft), sup(Iright))

    inside = (inf(Iθm) ≥ inf(Iallowed)) && (sup(Iθm) ≤ sup(Iallowed))

    if verbose
        println("θ_m ∈ ", Iθm, " ; allowed ∈ ", Iallowed, " ; inside? ", inside)
    end

    if !inside

        return true, :zero_roots, Iθm, interval(0), Iallowed
    end


    If3 = f₃(Iθm, Iκ, Iμ)
    tg = sign_tag(If3)

    if verbose
        println("f₃(θ_m) ∈ ", If3, " ⇒ ", tg)
    end

    if tg === :neg
        return true, :two_roots, Iθm, If3, Iallowed
    elseif tg === :pos
        return true, :zero_roots, Iθm, If3, Iallowed
    else
        return false, :unknown, Iθm, If3, Iallowed
    end
end

function θm_from_br(br, κ, μ; atol=1e-12, verbose=false)
    bestθ = NaN
    bestval = Inf
    found = false

    for x0 in br
        θ = find_zero(θ -> f₃θ(θ, κ, μ), x0; atol=atol)

        s2 = f₃θθ(θ, κ, μ)      # second derivative test
        if s2 > 0               # local minimum
            val = f₃(θ, κ, μ)   # height of the minimum
            if !found || val < bestval
                bestval = val
                bestθ = θ
                found = true
            end
        end

        verbose && println("x0=$x0  θ*=$θ  f₃θθ=$s2  f₃=$(f₃(θ, κ, μ))")
    end

    return found ? bestθ : NaN
end


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

    boxes = Tuple{Symbol,Float64,Float64,Float64,Float64,Bool,Bool,Symbol,Symbol}[]
    total = nx * ny
    seeded = 0
    validated = 0

    # iterate cell centers
    for i in 1:nx, j in 1:ny
        κc = κ_lo + (i - 0.5) * Δκ
        μc = μ_lo + (j - 0.5) * Δμ

        (θmin, θmax) = Intθ_3(κc, μc)

        br = brackets_signchanges(θ -> f₃θ(θ, κc, μc), θmin, θmax)

        if isempty(br)
            continue
        end
        seeded += 1 # Seeded means we found a root 

        θm = θm_from_br(br, κc, μc; atol=1e-12)
        # 2) refine the root 
        root = find_zero(θ -> f₃θ(θ, κc, μc), θm; atol=1e-12)

        Iκ = interval(κc, rκ; format=:midpoint)
        Iμ = interval(μc, rμ; format=:midpoint)

        θmax_refined = sup(Intθ_3(Iκ, Iμ)[2])
        θmin_refined = inf(Intθ_3(Iκ, Iμ)[1])

        ok_root, data = validate_root_sign_mono(
            f₃θ,
            f₃θθ,
            root,
            interval(θmin_refined, θmax_refined),
            Iκ,
            Iμ;
            ntiles=16,
            min_width=1e-12,
            min_radius=1e-12,
            verbose=false,
        )

        keep_invalid = true

        if ok_root

            clasif = classify_f3_rootcount(data.Iexist, Iκ, Iμ; verbose=false)

            if clasif[2] === :unknown
                push!(boxes, (:unknown, κc, μc, rκ, rμ, false, false, :unknown, :unknown))

            elseif clasif[2] === :zero_roots
                push!(boxes, (:zero_roots, κc, μc, rκ, rμ, false, false, :unknown, :unknown))
                validated += 1  

            else

                br = brackets_signchanges(θ -> f₃(θ, κc, μc), θmin, θmax)


                root1 = find_zero(θ -> f₃(θ, κc, μc), br[1]; atol=1e-12)
                root2 = find_zero(θ -> f₃(θ, κc, μc), br[2]; atol=1e-12)

                ### order roots by absolute value
                if abs(root2) < abs(root1)
                    (root1, root2) = (root2, root1)
                end

                ok_root1, data1 = validate_root_sign_mono(
                    f₃,
                    f₃θ,
                    root1,
                    interval(θmin_refined, θmax_refined),
                    Iκ,
                    Iμ;
                    ntiles=30,
                    min_width=1e-12,
                    min_radius=1e-12,
                    verbose=false,
                )

                ok_root2, data2 = validate_root_sign_mono(
                    f₃,
                    f₃θ,
                    root2,
                    interval(θmin_refined, θmax_refined),
                    Iκ,
                    Iμ;
                    ntiles=30,
                    min_width=1e-12,
                    min_radius=1e-12,
                    verbose=false,
                )

                 _, stability_label1 = classify_equilibrium(data1.Iexist, Iκ, Iμ)
                 _, stability_label2 = classify_equilibrium(data2.Iexist, Iκ, Iμ)




                push!(boxes, ( clasif[2], κc, μc, rκ, rμ, ok_root1, ok_root2, stability_label1, stability_label2))
                validated += 1
            end


        elseif keep_invalid

            push!(boxes, (:unknown, κc, μc, rκ, rμ, false, false, :unknown, :unknown))
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

nx, ny = 2 * 2464, 2 * 89
boxes, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


using CairoMakie

"""
    plot_regions_all_labels_from_boxes(boxes;
        label_type=:root,                 # :root or :stability
        root_index=1,                     # only used when label_type==:stability (1 or 2)
        κlims=(0.0, (0.5*pi)^2),
        μlims=(0.0, 0.0886765),
        resolution=2000,
        alpha=0.7,
        show_boundary=true,
        boundary_lw=2.0,
        title="",
        axis_box=true,
        unknown_box_lw=0.8)

Plots regions from a *regular box cover* by rasterizing the union of boxes per label.

NEW box structure expected:
    (rootcount::Symbol, κc::Float64, μc::Float64, rκ::Float64, rμ::Float64,
     ok1::Bool, ok2::Bool, stab1::Symbol, stab2::Symbol)

Behavior:
- label_type == :root:
    rasterizes regions for each rootcount label.
    (:unknown is NOT rasterized; its boxes are drawn at the end as black rectangles.)
- label_type == :stability:
    for boxes with rootcount == :zero_roots → rasterize as :zero_roots (light gray)
    else use stab1/stab2 depending on root_index.
    (:unknown stability is NOT rasterized; its boxes are drawn at the end as black rectangles.)

Returns: (fig, ax)
Call/signature and return are unchanged.
"""
function plot_regions_all_labels_from_boxes(boxes;
    label_type::Symbol = :root,
    root_index::Int = 1,
    κlims::Tuple{<:Real,<:Real} = (0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real} = (0.0, 0.0886765),
    resolution::Int = 2000,
    alpha::Real = 0.7,
    show_boundary::Bool = true,
    boundary_lw::Real = 2.0,
    title::AbstractString = "",
    axis_box::Bool = true,
    unknown_box_lw::Real = 0.8,
)

    κlo, κhi = κlims
    μlo, μhi = μlims
    @assert κhi > κlo
    @assert μhi > μlo
    @assert resolution ≥ 2
    @assert (label_type === :root) || (label_type === :stability)
    @assert (root_index == 1) || (root_index == 2)

    # grid endpoints (Makie likes explicit vectors for heatmap/contour)
    κgrid = collect(range(float(κlo), float(κhi); length=resolution))
    μgrid = collect(range(float(μlo), float(μhi); length=resolution))
    dκ = (κhi - κlo) / (resolution - 1)
    dμ = (μhi - μlo) / (resolution - 1)

    @inline clampi(i) = ifelse(i < 1, 1, ifelse(i > resolution, resolution, i))

    # ------------------- collect labels present (excluding :unknown; we draw unknown as boxes) -------------------
    labels = Symbol[]
    for bx in boxes
        rootcount, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        if label_type === :root
            rootcount === :unknown && continue
            push!(labels, rootcount)
        else
            if rootcount === :zero_roots
                push!(labels, :zero_roots)
            else
                stab = (root_index == 1) ? stab1 : stab2
                stab === :unknown && continue
                push!(labels, stab)
            end
        end
    end
    labels = sort!(unique(labels); by=String)

    # ------------------- allocate one mask per label -------------------
    masks = Dict{Symbol, Matrix{UInt8}}()
    for lbl in labels
        masks[lbl] = fill(UInt8(0), resolution, resolution)
    end

    # ------------------- rasterize each label's union -------------------
    for bx in boxes
        rootcount, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        # choose label for this box depending on mode
        # lbl::Symbol
        if label_type === :root
            rootcount === :unknown && continue
            lbl = rootcount
        else
            if rootcount === :zero_roots
                lbl = :zero_roots
            else
                stab = (root_index == 1) ? stab1 : stab2
                stab === :unknown && continue
                lbl = stab
                # optional gating (comment/uncomment as desired):
                # (root_index == 1 ? ok1 : ok2) || continue
            end
        end

        haskey(masks, lbl) || continue  # safety

        κa = κc - rκ; κb = κc + rκ
        μa = μc - rμ; μb = μc + rμ

        (κb < κlo || κa > κhi || μb < μlo || μa > μhi) && continue
        κa = max(κa, κlo); κb = min(κb, κhi)
        μa = max(μa, μlo); μb = min(μb, μhi)

        i1 = clampi(floor(Int, (κa - κlo) / dκ) + 1)
        i2 = clampi(ceil(Int,  (κb - κlo) / dκ) + 1)
        j1 = clampi(floor(Int, (μa - μlo) / dμ) + 1)
        j2 = clampi(ceil(Int,  (μb - μlo) / dμ) + 1)

        M = masks[lbl]
        @inbounds for i in i1:i2, j in j1:j2
            M[i, j] = UInt8(1)
        end
    end

    # ------------------- figure/axis -------------------
    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel = L"\mathbf{\kappa}",
        ylabel = L"\mathbf{\mu}",
        title  = title,
        ylabelsize = 24,
        xlabelsize = 24,
        xticklabelsize = 18,
        yticklabelsize = 18,
        topspinevisible = true,
        rightspinevisible = true,
        xgridvisible = false,
        ygridvisible = false,
        xtickalign = 1,
        ytickalign = 1,
        xticksmirrored = true,
        yticksmirrored = true,
        xticksize = 10,
        yticksize = 10,
        spinewidth = 2,
        xtickwidth = 2,
        ytickwidth = 2,
    )

    if axis_box
        ax.topspinevisible = true
        ax.rightspinevisible = true
    end

    # xlims!(ax, κlo, κhi)
    # ylims!(ax, μlo, μhi)

    # choose color function (you provide these)
    colorfn = (label_type === :root) ? rootcount_color : stability_color

    # ------------------- draw masks (rasterized unions) -------------------
    for lbl in labels
        M = masks[lbl]
        Z = Float32.(M)
        col = Makie.to_color(colorfn(lbl))
        α = float(alpha)

        heatmap!(ax, κgrid, μgrid, Z;
            colormap = [RGBAf(0,0,0,0), col],
            colorrange = (0f0, 1f0),
            interpolate = false,
            alpha = α,
        )

        if show_boundary
            contour!(ax, κgrid, μgrid, Z;
                levels = [0.5f0],
                color = col,
                linewidth = boundary_lw,
                labels = false,
            )
        end
    end

    # ------------------- draw UNKNOWN as actual boxes (black), on top -------------------
    # (This is what you asked: do NOT rasterize unknown pixels; draw the unknown boxes.)
    black = RGBf(0,0,0)

    for bx in boxes
        rootcount, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        is_unknown = if label_type === :root
            rootcount === :unknown
        else
            # only unknown stability where there is a root;
            # keep :zero_roots rasterized as gray.
            if rootcount === :zero_roots
                false
            else
                ((root_index == 1) ? stab1 : stab2) === :unknown
            end
        end
        is_unknown || continue

        κa = κc - rκ; κb = κc + rκ
        μa = μc - rμ; μb = μc + rμ

        (κb < κlo || κa > κhi || μb < μlo || μa > μhi) && continue
        κa = max(κa, κlo); κb = min(κb, κhi)
        μa = max(μa, μlo); μb = min(μb, μhi)

        # rectangle outline only (keeps file light compared to filled polys)
        lines!(ax,
            [κa, κb, κb, κa, κa],
            [μa, μa, μb, μb, μa];
            color = black,
            linewidth = unknown_box_lw
        )
    end

    return fig, ax
end


# --- Root count plot ---
fig1, ax1 = plot_regions_all_labels_from_boxes(boxes;
    label_type = :root,
    title = "Root-count",
)

save("figures/f3_root_regions.png", fig1)

# --- Stability plot (root 1) ---
fig2, ax2 = plot_regions_all_labels_from_boxes(boxes;
    label_type = :stability,
    root_index = 1,
    title = "Stability (root 1)",
)

save("figures/f3_stability_root1.png", fig2)

# --- Stability plot (root 2) ---
fig3, ax3 = plot_regions_all_labels_from_boxes(boxes;
    label_type = :stability,
    root_index = 2,
    title = "Stability (root 2)",
)

save("figures/f3_stability_root2.png", fig3)





# fig1, ax = plot_regions_all_labels_from_boxes(
#     boxes;
#     label_type=:root,
#     κlims=(0.0, (0.5 * pi)^2),
#     μlims=(0.0, 0.0886765),
#     resolution=2000,
#     alpha=0.7,
#     show_boundary=true,
#     boundary_lw=2.0,
#     title="Root-count Regions",
# )


# CairoMakie.save("figures/f3_root_regions.png", fig1)
