



κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 2000, 2000
boxes1 = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


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


κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 2000, 2000
boxes2 = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


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


κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 2000, 2000
boxes3 = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


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


