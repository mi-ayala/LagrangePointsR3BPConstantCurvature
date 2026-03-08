

#### Negative and positive curvature together
include("numerics_f5.jl")
include("numerics_f5neg.jl")


κ_lo, κ_hi = -1.0, 1.0
μ_lo, μ_hi = 0.0, 0.0886765

append!(boxes5, boxes5neg)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    κlims=(κ_lo, κ_hi),
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f5_stability_root1_negative_numerics_together.png", fig2)

