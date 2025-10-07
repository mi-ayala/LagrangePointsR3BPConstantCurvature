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

function f₁θ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r  = sqrt(κ)
    a  = θ + r*q₁(κ, μ)
    b  = θ - r*q₂(κ, μ)
    Γv = Γ(κ, μ)

    d_main = cos(twice(θ))

    sa = sin(a);  ca = cos(a)
    sb = sin(b);  cb = cos(b)

    # d/dθ [ 1/sin(a)^2 - μ/sin(b)^2 ]
    #   = -2 cos(a)/sin(a)^3 + 2 μ cos(b)/sin(b)^3
    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)
    d_term  = (-twice(ca) * inv_sa3) + (twice(μ) * cb * inv_sb3)

    return d_main - Γv * d_term
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
@show f₁θ(interval(0.1), κ, μ)



using ForwardDiff

κ = 0.5
μ = 0.02
θ = 0.7

# compare analytic vs automatic derivative
ad = ForwardDiff.derivative(θ -> f₁(θ, κ, μ), θ)
an = f₁θ(θ, κ, μ)

println("analytic  f₁θ = ", an)
println("ForwardDiff = ", ad)
println("difference  = ", ad - an)



# bounds
μT = 0.0886765
κ_lo, κ_hi = 0.0, (0.5*pi)^2       # (0, π/2)
μ_lo, μ_hi = 0.0, μT           # (0, μT)

# make an interior grid (avoid exact endpoints to stay away from singular boundaries)
interior_linspace(lo, hi, n; margin=0.02) = begin
    a = lo + margin*(hi - lo)
    b = hi - margin*(hi - lo)
    range(a, b; length=n)
end

nk, nm = 15, 15
κs = collect(interior_linspace(κ_lo, κ_hi, nk))  # κ-grid
μs = collect(interior_linspace(μ_lo, μ_hi, nm))  # μ-grid

# Example: make all grid pairs as a vector of tuples (κ, μ)
grid = [(κ, μ) for κ in κs, μ in μs]  # nk×nm matrix of pairs

# Optionally flatten (e.g. for map):
pairs = vec(grid)  # length = nk*nm



using Roots


function brackets(f, a::Real, b::Real; n::Integer=600)
    xs = collect(range(a, b; length=n+1))
    br = Tuple{Float64,Float64}[]
    fa = f(xs[1])
    for i in 1:n
        x1, x2 = xs[i], xs[i+1]
        fx = f(x2)
        if isfinite(fa) && isfinite(fx) && sign(fa) != sign(fx)
            push!(br, (float(x1), float(x2)))
        end
        fa = fx
    end
    br
end

function roots_brent_f₁(κ, μ; n=600, digits=12, shrink=0.3)
    θmin = -sqrt(κ)*q₁(κ, μ) + shrink
    θmax =  sqrt(κ)*q₂(κ, μ) - shrink
    f(θ) = f₁(θ, κ, μ)

    rts = Float64[]
    for (a,b) in brackets(f, θmin, θmax; n=n)
        try
            r = find_zero(f, (a,b), Brent())
            isfinite(r) && push!(rts, r)
        catch
            # skip failed brackets (often near poles)
        end
    end

    # normalize by rounding, then sort + dedup
    rts = round.(rts; digits=digits)
    sort!(rts)
    # deduplicate adjacent equal values after rounding
    if !isempty(rts)
        dedup = Float64[rts[1]]
        for x in rts[2:end]
            x == dedup[end] || push!(dedup, x)
        end
        rts = dedup
    end
    return rts
end




κ, μ =2.0, 0.07
rts = roots_brent_f₁(κ, μ)  # numeric roots in (θmin, θmax)


Intθ(κ, μ)

using Plots
Plots.plot(θ -> f₁(θ, κ, μ), 0.25, 1)

for i in range(Intθ(κ, μ)[1],Intθ(κ, μ)[2]; length=400)
    println("f₁($i) = ", f₁(i, κ, μ))
end