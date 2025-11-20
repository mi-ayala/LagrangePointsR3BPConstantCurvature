include("functions.jl")


μ= 0.05

# f= (u, γ) -> f₂(u, γ, μ)
# Dᵤf=(u, γ) -> f₂θ(u, γ, μ)
# ∂γf=(u, γ) -> f₂μ(u, κ, γ) 


f = (u, γ) -> [f₂(u, κ, μ); f₂θ(u, κ, μ)]



import LinearAlgebra: ⋅

F(x, v, w) = [f(x[1:2], x[3]) ; (x - w) ⋅ v]

DF(x, v) = [Dᵤf(x[1:2], x[3]) ∂γf(x[1:2], x[3]) ; transpose(v)]


using RadiiPolynomial
import LinearAlgebra: nullspace

# initial point on the branch of equilibria

γ_init = 2.0

u_init = [1.1, 0.5]
u_init, success = newton(u -> (f(u, γ_init), Dᵤf(u, γ_init)), u_init)

# next point on the branch of equilibria

x_init = [u_init ; γ_init]

v = vec(nullspace([Dᵤf(x_init[1:2], x_init[3]) ∂γf(x_init[1:2], x_init[3])])) # initial tangent vector

δ = 5e-2 # step size

w = x_init + δ * v # predictor

x_final, success = newton(x -> (F(x, v, w), DF(x, v)), w)
