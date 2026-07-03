using RC3BP
using CairoMakie
using LaTeXStrings

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
save(joinpath(@__DIR__, "..", "figures", "numerics", "f1_root_regions_numerics.png"), fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes1;
    label_type=:stability,
    root_index=1,
    ms=4
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "f1_stability_numerics.png"), fig2)


### Plots f2
fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes2;
    label_type=:root,
    ms=4,
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "f2_root_regions_numerics.png"), fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes2;
    label_type=:stability,
    root_index=1,
    ms=4
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "f2_stability_root1_numerics.png"), fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes2;
    label_type=:stability,
    root_index=2,
    ms=4
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "f2_stability_root2_numerics.png"), fig3)


### Plots f3

fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes3;
    label_type=:root,
    ms=4,
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "f3_root_regions_numerics.png"), fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes3;
    label_type=:stability,
    root_index=1,
    ms=4
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "f3_stability_root1_numerics.png"), fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes3;
    label_type=:stability,
    root_index=2,
    ms=4
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "f3_stability_root2_numerics.png"), fig3)


### Plots f4

fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes4;
    label_type=:root,
    ms=4
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "f4_root_regions_numerics.png"), fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes4;
    label_type=:stability,
    root_index=1,
    ms=4
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "f4_stability_numerics.png"), fig2)


#### Function to count roots through all the functions.
boxes_total = combine_rootcount_boxes([boxes1, boxes2, boxes3, boxes4])



fig, ax = plot_regions_all_labels_from_total_boxes(boxes_total;
    ms=4,
    # title="Total number of roots",
    κlims=(-1.0, (0.5 * pi)^2),
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "total_root_regions.png"), fig)



#### Bifurcation diagram
fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa_collinear(boxes, 0.01;
    use=:stability,
    ms=5,
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "collinear_bif_plot_mu01.png"), fig3)



fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa_collinear(boxes, 0.06;
    use=:stability,
    ms=5,
)
save(joinpath(@__DIR__, "..", "figures", "numerics", "collinear_bif_plot_mu06.png"), fig3)
