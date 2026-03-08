
include("functions_plots.jl")
include("numerics_f5.jl")
include("numerics_f5neg.jl")



# Plot
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



#### Together with negative κ
append!(boxes5, boxes5neg)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    κlims=(-(0.5 * pi)^2, (0.5 * pi)^2),
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f5_stability_root1_negative_numerics_together.png", fig2)
