include("numerics_f1.jl")
include("numerics_f2.jl")
include("numerics_f3.jl")
include("numerics_f4.jl")


using CairoMakie
function waterfall_regions_from_boxes_clean(boxes_list;
    names=nothing,
    label_type::Symbol = :root,
    root_index::Int = 1,
    κlims::Tuple{<:Real,<:Real} = (0.0, (0.5*pi)^2),
    μlims::Tuple{<:Real,<:Real} = (0.0, 0.0886765),
    ms::Real = 3.0,
    alpha::Real = 0.05,          # translucent
    yspacing::Real = 1.0,
    title::AbstractString = "Waterfall regions",
)

    κlo, κhi = κlims
    μlo, μhi = μlims
    nsets = length(boxes_list)

    if names === nothing
        names = ["f$i" for i in 1:nsets]
    end

    colorfn = (label_type === :root) ? rootcount_color : stability_color

    fig = Figure()
    ax = Axis3(fig[1,1],
        title=title,
        xlabel=L"\kappa",
        ylabel="",
        zlabel=L"\mu",
    )

    xlims!(ax, κlo, κhi)
    zlims!(ax, μlo, μhi)
    ylims!(ax, 0, (nsets-1)*yspacing)

    for (i, boxes) in enumerate(boxes_list)
        y0 = (i-1)*yspacing

        pts_by_lbl = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()

        for bx in boxes
            root_lbl, κc, μc, θ1, θ2, stab1, stab2 = bx

            (κc < κlo || κc > κhi || μc < μlo || μc > μhi) && continue

            if label_type === :root
                # skip zero-root slices completely
                root_lbl === :zero_roots && continue
                lbl = root_lbl
            else
                # stability: skip if this root does not exist
                θ = (root_index == 1) ? θ1 : θ2
                isfinite(θ) || continue
                lbl = (root_index == 1) ? stab1 : stab2
            end

            xy = get!(pts_by_lbl, lbl) do
                (Float64[], Float64[])
            end
            push!(xy[1], κc)
            push!(xy[2], μc)
        end

        for (lbl, (κvals, μvals)) in pts_by_lbl
            isempty(κvals) && continue

            col = Makie.to_color(colorfn(lbl))

            scatter!(ax,
                κvals,
                fill(y0, length(κvals)),
                μvals;
                markersize=ms,
                color=col,
                alpha=alpha,
                transparency=true,
            )
        end
    end

    hidespines!(ax)
    # ax.gridvisible = false

    return fig, ax
end

figW, axW = waterfall_regions_from_boxes_clean(
    [boxes1, boxes2, boxes3, boxes4];
    names=["f1","f2","f3","f4"],
    label_type=:root,
    title="",
)

save("figures/numerics/waterfall_root_regions.png", figW)