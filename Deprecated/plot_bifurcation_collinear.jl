

include("numerics_f1.jl")
include("numerics_f2.jl")
include("numerics_f3.jl")
include("numerics_f4.jl")


boxes = Any[]
append!(boxes, boxes1)
append!(boxes, boxes2)
append!(boxes, boxes3)
append!(boxes, boxes4)



fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa(boxes, 0.01;
    use=:stability,
    ms=5,
)

save("figures/numerics/collinear_bif_plot_mu01.png", fig3)



fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa(boxes, 0.06;
    ms=5,
)

save("figures/numerics/collinear_bif_plot_mu06.png", fig3)

