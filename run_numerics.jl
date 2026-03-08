



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