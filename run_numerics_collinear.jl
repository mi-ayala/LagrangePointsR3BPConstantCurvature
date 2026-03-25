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
    label_type = :root,
    ms = 4,
)

save("figures/numerics/f2_root_regions_numerics.png", fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes2;
    label_type = :stability,
    root_index = 1,
        ms = 4
)

save("figures/numerics/f2_stability_root1_numerics.png", fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes2;
    label_type = :stability,
    root_index = 2,
        ms = 4
)

save("figures/numerics/f2_stability_root2_numerics.png", fig3)


### Plots f3

fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes3;
    label_type = :root,
    ms = 4,
)

save("figures/numerics/f3_root_regions_numerics.png", fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes3;
    label_type = :stability,
    root_index = 1,
        ms = 4
)

save("figures/numerics/f3_stability_root1_numerics.png", fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes3;
    label_type = :stability,
    root_index = 2,
        ms = 4
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
        lbl === :zero_roots  ? 0 :
        lbl === :one_root    ? 1 :
        lbl === :two_roots   ? 2 :
        lbl === :three_roots ? 3 :
        lbl === :four_roots  ? 4 :
        lbl === :five_roots  ? 5 :
        lbl === :six_roots   ? 6 :
        lbl === :unknown     ? -1 :
                               error("Unknown root-count label: $lbl")
    )

    # map counts -> labels
    countlabel(n::Int) = (
        n == 0 ? :zero_roots  :
        n == 1 ? :one_root    :
        n == 2 ? :two_roots   :
        n == 3 ? :three_roots :
        n == 4 ? :four_roots  :
        n == 5 ? :five_roots  :
        n == 6 ? :six_roots   :
                 error("Total root count $n not supported")
    )

    # build a dictionary for each boxes array: (κ,μ) -> label
    dicts = Vector{Dict{Tuple{Float64,Float64},Symbol}}()
    for boxes in boxes_list
        d = Dict{Tuple{Float64,Float64},Symbol}()
        for bx in boxes
            lbl = bx[1]
            κc  = float(bx[2])
            μc  = float(bx[3])
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

fig, ax = plot_regions_all_labels_from_total_boxes(boxes_total;
    ms=4,
    # title="Total number of roots",
    κlims=(-1.0, (0.5 * pi)^2),
)

save("figures/numerics/total_root_regions.png", fig)



#### Bifurcation diagram

fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa(boxes, 0.01;
    use=:stability,
    ms=5,
)

save("figures/numerics/collinear_bif_plot_mu01.png", fig3)



fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa(boxes, 0.06;
    ms=5,
)

save("figures/numerics/collinear_bif_plot_mu06.png", fig3)
