using IntervalArithmetic
using ForwardDiff

# include("functions.jl")
include("functions_f4.jl")

### Testing interval IntervalArithmetic
κ = interval(0.5)
μ = interval(0.5)
# @show q₁(κ, μ)
# @show q₂(κ, μ)
# @show Γ(κ, μ)
# @show f₁(interval(0.1), κ, μ)
# @show Intθ_1(κ, μ)
# @show f₁θ(interval(0.1), κ, μ)
# @show f₂(interval(0.1), κ, μ)
# @show Intθ_2(κ, μ)
# @show f₃(interval(0.1), κ, μ)
# @show Intθ_3(κ, μ)
# @show f₄(interval(0.1), κ, μ)
# @show Intθ_4(κ, μ)
# @show Δ(interval(0.1), κ, μ)
# @show c₂(interval(0.1), κ, μ)
# @show f₂θ(interval(0.1), κ, μ)
# @show f₂κ(interval(0.1), κ, μ)
@show f₄(interval(0.1), κ, μ)
@show f₄θ(interval(0.1), κ, μ)



 # compare analytic vs automatic derivative
κ = 0.5
μ = 0.2
θ = 0.7
# ad = ForwardDiff.derivative(θ -> f₂(θ, κ, μ), θ)
# an = f₂θ(θ, κ, μ)
# println("analytic  f₄θ = ", an)
# println("ForwardDiff = ", ad)
# println("difference  = ", ad - an)

ad = ForwardDiff.derivative(θ -> f₄(θ, κ, μ), θ)
an = f₄θ(θ, κ, μ)
println("analytic  f₄θ = ", an)
println("ForwardDiff = ", ad)
println("difference  = ", ad - an)


# ### kappa derivative
# ad = ForwardDiff.derivative(κ -> f₂(0.7, κ, μ), κ)
# an = f₂κ(0.7, κ, μ)
# println("analytic  f₂κ = ", an)
# println("ForwardDiff = ", ad)
# println("difference  = ", ad - an)


# ### mu derivative
# ad = ForwardDiff.derivative(μ -> f₂(0.7, κ, μ), μ)
# an = f₂μ(0.7, κ, μ)
# println("analytic  f₂μ = ", an)
# println("ForwardDiff = ", ad)
# println("difference  = ", ad - an)  

# # compare analytic vs automatic derivative
# κ = 0.5
# μ = 0.02
# θ = 0.7
# ad = ForwardDiff.derivative(θ -> f₁(θ, κ, μ), θ)
# an = f₁θ(θ, κ, μ)

# println("analytic  f₁θ = ", an)
# println("ForwardDiff = ", ad)
# println("difference  = ", ad - an)


# ### A single root Example 
# κ, μ = 0.5, 0.02
# ### Estimate by sign changes
# for i in range(Intθ(κ, μ)[1],Intθ(κ, μ)[2]; length=400)
#     println("f₁($i) = ", f₁(i, κ, μ))
# end


# br = brackets_f₁(κ, μ; n=800)
# f(θ) = f₁(θ, κ, μ)
# rts = refine_brackets(f, br)
# println(rts)
# find_zero(θ -> f₁(θ, κ, μ), br[1])


# br = brackets_f₂(κ, μ; n=800)

# rts = refine_brackets(f, br)
# println(rts)
# find_zero(θ -> f₁(θ, κ, μ), br[1])


# ### Testing plotting function
# boxes = [
#     (true, 0.02, 0.02, 1e-3, 2e-3),
#     # skipped
# ]
# plot_validated_boxes(boxes; color=:orange, alpha=0.5, title="Validated regions in (κ, μ)")

# # ------------------------------------------------------------
# # Quick test block: call all derivatives with Interval inputs
# # ------------------------------------------------------------

# using IntervalArithmetic

# # Example parameter values (replace by whatever you want)
# κ  = interval(0.1)
# μ  = interval(0.02)
# θ  = interval(0.3)

# println("------ Testing f₂ and all derivatives on intervals ------")

# @show f₂(θ, κ, μ)

# @show f₂θ(θ, κ, μ)
# @show f₂κ(θ, κ, μ)
# @show f₂μ(θ, κ, μ)

# @show f₂θθ(θ, κ, μ)
# @show f₂θκ(θ, κ, μ)
# @show f₂θμ(θ, κ, μ)

# println("------ Testing q₁, q₂ and their derivatives ------")

# @show q₁(κ, μ)
# @show q₂(κ, μ)

# @show q₁κ(κ, μ)
# @show q₁μ(κ, μ)

# println("------ Testing Γ and its derivatives ------")

# @show Γ(κ, μ)
# @show Γκ(κ, μ)
# @show Γμ(κ, μ)

# using ForwardDiff
# using StaticArrays   # for a tiny F(u) definition

# # Pick a test point away from singularities (adjust if needed)
# θ0 = 0.35
# κ0 = 0.7
# μ0 = 0.01

# println("=== First derivatives of f₂ ===")

# ad_f2θ = ForwardDiff.derivative(t -> f₂(t, κ0, μ0), θ0)
# an_f2θ = f₂θ(θ0, κ0, μ0)
# @show an_f2θ ad_f2θ an_f2θ - ad_f2θ

# ad_f2κ = ForwardDiff.derivative(k -> f₂(θ0, k, μ0), κ0)
# an_f2κ = f₂κ(θ0, κ0, μ0)
# @show an_f2κ ad_f2κ an_f2κ - ad_f2κ

# ad_f2μ = ForwardDiff.derivative(m -> f₂(θ0, κ0, m), μ0)
# an_f2μ = f₂μ(θ0, κ0, μ0)
# @show an_f2μ ad_f2μ an_f2μ - ad_f2μ

# println("\n=== Second / mixed derivatives via f₂θ ===")

# ad_f2θθ = ForwardDiff.derivative(t -> f₂θ(t, κ0, μ0), θ0)
# an_f2θθ = f₂θθ(θ0, κ0, μ0)
# @show an_f2θθ ad_f2θθ an_f2θθ - ad_f2θθ

# ad_f2θκ = ForwardDiff.derivative(k -> f₂θ(θ0, k, μ0), κ0)
# an_f2θκ = f₂θκ(θ0, κ0, μ0)
# @show an_f2θκ ad_f2θκ an_f2θκ - ad_f2θκ

# ad_f2θμ = ForwardDiff.derivative(m -> f₂θ(θ0, κ0, m), μ0)
# an_f2θμ = f₂θμ(θ0, κ0, μ0)
# @show an_f2θμ ad_f2θμ an_f2θμ - ad_f2θμ

# println("\n=== Jacobian J_f₂ vs ForwardDiff.jacobian ===")

# # Analytic Jacobian
# J_an = J_f₂(θ0, κ0, μ0)
# @show J_an

# # AD Jacobian of F(θ, κ, μ) = [f₂; f₂θ]
# F(u) = @SVector [ f₂(u[1], u[2], u[3]),
#                   f₂θ(u[1], u[2], u[3]) ]

# u0   = SVector(θ0, κ0, μ0)
# J_ad = ForwardDiff.jacobian(F, u0)
# @show J_ad

# println("\nDifference J_an - J_ad:")
# @show J_an .- J_ad