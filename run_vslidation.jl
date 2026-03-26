using RC3BP

include("validation_f1.jl")
include("validation_f2.jl")
include("validation_f3.jl")
include("validation_f4.jl")

boxes = Any[]
append!(boxes, boxes1_validated)
append!(boxes, boxes2_validated)
append!(boxes, boxes3_validated)
append!(boxes, boxes4_validated)


function plot_unknown_boxes(boxes;
    label_type::Symbol=:root,
    root_index::Int=1,
    κlims::Tuple{<:Real,<:Real}=(0.0, (0.5 * pi)^2),
    μlims::Tuple{<:Real,<:Real}=(0.0, 0.0886765),
    resolution::Int=2000,
    alpha::Real=1.0,
    title::AbstractString="",
    axis_box::Bool=true,

    # --- pretty analytical region overlay ---
    overlay_rect=nothing,                    # (κa, κb, μa, μb)
    overlay_nx::Int=0,
    overlay_ny::Int=0,
    overlay_color=RGBAf(0.97, 0.98, 1.00, 1),
    overlay_line_color=RGBAf(1, 1, 1, 1),
    overlay_line_alpha::Real=0.95,
    overlay_line_lw::Real=0.9,
    overlay_bg_alpha::Real=0.30,
)
    κlo, κhi = κlims
    μlo, μhi = μlims

    @assert κhi > κlo
    @assert μhi > μlo
    @assert resolution ≥ 2
    @assert (label_type === :root) || (label_type === :stability)
    @assert root_index == 1 || root_index == 2

    κgrid = collect(range(float(κlo), float(κhi); length=resolution))
    μgrid = collect(range(float(μlo), float(μhi); length=resolution))
    dκ = (κhi - κlo) / (resolution - 1)
    dμ = (μhi - μlo) / (resolution - 1)

    @inline clampi(i) = ifelse(i < 1, 1, ifelse(i > resolution, resolution, i))

    # --- unknown mask only ---
    M = fill(UInt8(0), resolution, resolution)

    # --- rasterize only unknown boxes ---
    for bx in boxes
        root_lbl, κc, μc, rκ, rμ, ok1, ok2, stab1, stab2 = bx

        lbl = if label_type === :root
            root_lbl
        else
            root_lbl === :zero_roots ? :zero_roots :
            (root_index == 1 ? stab1 : stab2)
        end

        lbl === :unknown || continue

        κa = max(κc - rκ, κlo)
        κb = min(κc + rκ, κhi)
        μa = max(μc - rμ, μlo)
        μb = min(μc + rμ, μhi)

        i1 = clampi(floor(Int, (κa - κlo) / dκ) + 1)
        i2 = clampi(ceil(Int, (κb - κlo) / dκ) + 1)
        j1 = clampi(floor(Int, (μa - μlo) / dμ) + 1)
        j2 = clampi(ceil(Int, (μb - μlo) / dμ) + 1)

        @inbounds for i in i1:i2, j in j1:j2
            M[i, j] = UInt8(1)
        end
    end

    # --- figure ---
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

    xlims!(ax, κlo, κhi)
    ylims!(ax, μlo, μhi)

    Z = Float32.(M)
    heatmap!(ax, κgrid, μgrid, Z;
        colormap=[RGBAf(0, 0, 0, 0), RGBf(0, 0, 0)],
        colorrange=(0f0, 1f0),
        interpolate=false,
        alpha=1.0,
    )

    # --- analytical region overlay ---
    if overlay_rect !== nothing
        κa, κb, μa, μb = overlay_rect
        κa = max(float(κa), κlo)
        κb = min(float(κb), κhi)
        μa = max(float(μa), μlo)
        μb = min(float(μb), μhi)

        if κb > κa && μb > μa
            poly!(ax,
                Point2f[(κa, μa), (κb, μa), (κb, μb), (κa, μb)];
                color=Makie.to_color(overlay_color),
                alpha=overlay_bg_alpha,
                strokewidth=0,
            )

            if overlay_nx > 0
                Δκg = (κb - κa) / overlay_nx
                for m in 1:(overlay_nx-1)
                    κx = κa + m * Δκg
                    lines!(ax, [κx, κx], [μa, μb];
                        color=overlay_line_color,
                        alpha=overlay_line_alpha,
                        linewidth=overlay_line_lw,
                    )
                end
            end

            if overlay_ny > 0
                Δμg = (μb - μa) / overlay_ny
                for m in 1:(overlay_ny-1)
                    μy = μa + m * Δμg
                    lines!(ax, [κa, κb], [μy, μy];
                        color=overlay_line_color,
                        alpha=overlay_line_alpha,
                        linewidth=overlay_line_lw,
                    )
                end
            end

            lines!(ax, [κa, κa], [μa, μb];
                color=overlay_line_color,
                linewidth=1.4,
                overdraw=true,
            )
            lines!(ax, [κb, κb], [μa, μb];
                color=overlay_line_color,
                linewidth=1.4,
                overdraw=true,
            )
            lines!(ax, [κa, κb], [μa, μa];
                color=overlay_line_color,
                linewidth=1.4,
                overdraw=true,
            )
            lines!(ax, [κa, κb], [μb, μb];
                color=overlay_line_color,
                linewidth=1.4,
                overdraw=true,
            )
        end
    end

    return fig, ax
end


fig, ax = plot_unknown_boxes(boxes;
    label_type=:root,
    title="",
    overlay_rect=(0.0, (pi^2 / 64), 0.0, 0.0886765),
    overlay_nx=10,
    overlay_ny=10,
)