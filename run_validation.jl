using RC3BP

include("functions_plots.jl")
include("validation_f1.jl")
include("validation_f2.jl")
include("validation_f3.jl")
include("validation_f4.jl")


boxes = Any[]
append!(boxes, boxes1_validated)
append!(boxes, boxes2_validated)
append!(boxes, boxes3_validated)
append!(boxes, boxes4_validated)


fig, ax = plot_unknown_boxes(boxes;
    label_type=:root,
    title="",
    overlay_rect=(0.0, (pi^2 / 64), 0.0, 0.0886765),
    overlay_nx=10,
    overlay_ny=10,
)

save("figures/validated/proof_regions.png", fig)

