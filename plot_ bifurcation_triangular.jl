include("numerics_f5.jl")
include("numerics_f5neg.jl")

boxes = Any[]
append!(boxes, boxes5)
append!(boxes, boxes5neg)

fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa(boxes, 0.01;
    use=:stability,
    ms=4,
    κlims=(0.0, 1.0)
)

save("figures/numerics/triangular_bif_plot_mu01.png", fig3)

fig3, ax3, μstar = plot_mu_slice_theta_over_sqrtkappa(boxes, 0.06;
    use=:stability,
    ms=4,
    κlims=(0.0, 1.0)
)

save("figures/numerics/triangular_bif_plot_mu06.png", fig3)

