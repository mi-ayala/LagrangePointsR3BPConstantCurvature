using IntervalArithmetic

twice(x) = x + x
sqr(x)   = x*x

function q₁(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    y = μ * sin(twice(r))                   
    x = one(μ) + μ * cos(twice(r))          
    return (one(r) / twice(r)) * atan(y, x) 
end

function q₂(κ, μ)
    κ, μ = promote(κ, μ)
    return one(κ) - q₁(κ, μ)
end

function Γ(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    num = sqr(sin(r)) * sin(r) * cos(r)     
    den = sqrt( μ*μ + twice( μ * cos(twice(r)) ) + one(μ) )  
    return num / den
end

function f₁(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r  = sqrt(κ)
    a  = θ + r*q₁(κ, μ)
    b  = θ - r*q₂(κ, μ)
    term = inv(sqr(sin(a))) - μ * inv(sqr(sin(b)))
    return sin(θ)*cos(θ) - Γ(κ, μ) * term
end

function Intθ(κ, μ)
    return (-sqrt(κ)*q₁(κ, μ), sqrt(κ)*q₂(κ, μ))
end


κ=interval(0.5)
μ=interval(0.5) 

@show q₁(κ, μ)
@show q₂(κ, μ)
@show Γ(κ, μ)
@show f₁(interval(0.1), κ, μ)
@show Intθ(κ, μ)