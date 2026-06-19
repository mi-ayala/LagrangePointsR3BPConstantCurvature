

include("numerics_f5.jl")
include("numerics_f5neg.jl")
append!(boxes5, boxes5neg)


include("functions_plots.jl")
# # Plot
fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    label_type=:root,
    ms=4,
)

save("figures/numerics/f5_root_regions_numerics2.png", fig1)


fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f5_stability_root1_numerics2.png", fig2)


fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes5;
    label_type=:stability,
    root_index=2,
    ms=4
)

save("figures/numerics/f5_stability_root2_numerics2.png", fig3)




# fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa_triangular(boxes5, 0.01;
#     use=:stability,
#     ms=5,
#     κlims=(-1.0, 1.0)
# )

# save("figures/numerics/triangular_bif_plot_mu01.png", fig3)

# fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa_triangular(boxes5, 0.06;
#     use=:stability,
#     ms=5,
#     κlims=(-1.0, 1.0)
# )

# save("figures/numerics/triangular_bif_plot_mu06.png", fig3)
