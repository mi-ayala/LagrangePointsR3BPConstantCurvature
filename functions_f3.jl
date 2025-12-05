function f₃(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) + μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) + Γ(κ, μ) * term
end

function Intθ_3(κ, μ)
    return (sqrt(κ) * q₂(κ, μ) - pi, -sqrt(κ) * q₁(κ, μ))
end

function brackets_f₃(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_3(κ, μ)[1]
    θmax = Intθ_3(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₃(θ, κ, μ), a, b; n=n)
end