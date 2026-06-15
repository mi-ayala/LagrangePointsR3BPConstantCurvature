using RC3BP
using JLD2

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

using JLD2

function proof_box(boxes)
    κmin = Inf
    κmax = -Inf
    μmin = Inf
    μmax = -Inf

    for (_, κc, μc, rκ, rμ, ok1, ok2, _, _) in boxes
        (ok1 || ok2) || continue

        κmin = min(κmin, κc - rκ)
        κmax = max(κmax, κc + rκ)
        μmin = min(μmin, μc - rμ)
        μmax = max(μmax, μc + rμ)
    end

    return (κmin, κmax, μmin, μmax)
end

proof_box1 = proof_box(boxes1_validated)
proof_box2 = proof_box(boxes2_validated)
proof_box3 = proof_box(boxes3_validated)
proof_box4 = proof_box(boxes4_validated)

println("boxes1_validated => ", proof_box1)
println("boxes2_validated => ", proof_box2)
println("boxes3_validated => ", proof_box3)
println("boxes4_validated => ", proof_box4)

@save "proof_boxes.jld2" proof_box1 proof_box2 proof_box3 proof_box4