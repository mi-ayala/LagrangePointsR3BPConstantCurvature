using IntervalArithmetic
using ForwardDiff

include("functions.jl")

### Testing interval IntervalArithmetic
κ = interval(0.5)
μ = interval(0.5)
@show q₁(κ, μ)
@show q₂(κ, μ)
@show Γ(κ, μ)
@show f₁(interval(0.1), κ, μ)
@show Intθ_1(κ, μ)
@show f₁θ(interval(0.1), κ, μ)
@show f₂(interval(0.1), κ, μ)
@show Intθ_2(κ, μ)
@show f₃(interval(0.1), κ, μ)
@show Intθ_3(κ, μ)
@show f₄(interval(0.1), κ, μ)
@show Intθ_4(κ, μ)
@show Δ(interval(0.1), κ, μ)
@show c₂(interval(0.1), κ, μ)
@show f₂θ(interval(0.1), κ, μ)
@show f₂κ(interval(0.1), κ, μ)
@show f₂μ(interval(0.1), κ, μ)



# compare analytic vs automatic derivative
κ = 0.5
μ = 0.2
θ = 0.7
ad = ForwardDiff.derivative(θ -> f₂(θ, κ, μ), θ)
an = f₂θ(θ, κ, μ)
println("analytic  f₄θ = ", an)
println("ForwardDiff = ", ad)
println("difference  = ", ad - an)


### kappa derivative
ad = ForwardDiff.derivative(κ -> f₂(0.7, κ, μ), κ)
an = f₂κ(0.7, κ, μ)
println("analytic  f₂κ = ", an)
println("ForwardDiff = ", ad)
println("difference  = ", ad - an)


### mu derivative
ad = ForwardDiff.derivative(μ -> f₂(0.7, κ, μ), μ)
an = f₂μ(0.7, κ, μ)
println("analytic  f₂μ = ", an)
println("ForwardDiff = ", ad)
println("difference  = ", ad - an)  

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