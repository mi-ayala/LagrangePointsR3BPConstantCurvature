module RC3BP

using CairoMakie
using IntervalArithmetic
using LaTeXStrings


include("functions_helpers.jl")
include("functions_numerics.jl")
include("functions_validation.jl")
include("functions_f1.jl")
include("functions_f1neg.jl")
include("functions_f2.jl")
include("functions_f2neg.jl")
include("functions_f3.jl")
include("functions_f3neg.jl")
include("functions_f4.jl")
include("functions_f4neg.jl")
include("functions_f5.jl")
include("functions_f5neg.jl")
include("functions_plots.jl")

export twice, csc, powabs3, sqr, three, inv3, cbrt_t
export q₁, q₂, Γ, λ₁, λ₂, detλ, c₂, Δ, q₁μ, q₁κ, Γμ, Γκ, Λ, Λκ, Λμ
export brackets_signchanges, classify_equilibrium_f, classify_center_equilibrium_num
export ispos, isneg, nonzero, eval_subdiv, certify_sign_interval_centered
export find_sign_change_centered, find_sign_change_centered_refine
export validate_root_sign_mono, sign_tag, classify_equilibrium
export complements_nonzero, finalize_root_label, complements_nonzero2
export classify_center_equilibrium
export f₁, f₁θ, Intθ_1, θ1_max, brackets_f₁
export f₁neg, f₁negθ
export f₂, f₂θ, f₂θθ, f₂κ, f₂μ, Intθ_2, brackets_f₂
export f₂neg, f₂negθ, f₂negθθ, f₂negκ, f₂negμ
export f₃, f₃θ, f₃θθ, Intθ_3, brackets_f₃
export f₃neg, f₃negθ, f₃negθθ
export f₄, f₄θ, Intθ_4, θ4_min, θ4_max, brackets_f₄
export f₄neg, f₄negθ
export sinc_t, f₅, f₅θ, Intθ_5
export q₁neg, q₂neg, Γneg, Λneg, f₅neg, classify_center_equilibrium_num_neg, Intθ_5neg
export rootcount_color, stability_color
export plot_regions_all_labels_from_boxes
export plot_regions_all_labels_from_boxes_numerics
export plot_regions_all_labels_from_total_boxes
export plot_mu_slice_theta_over_sqrtkappa

end # module RC3BP
