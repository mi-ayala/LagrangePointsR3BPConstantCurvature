using IntervalArithmetic
using ForwardDiff
using Roots
using Plots
using PlotlyJS
plotlyjs()


twice(x) = x + x
sqr(x) = x * x

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
    den = sqrt(μ * μ + twice(μ * cos(twice(r))) + one(μ))
    return num / den
end

### Fucntions and interval
function f₁(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) - μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) - Γ(κ, μ) * term
end

function Intθ_1(κ, μ)
    return (-sqrt(κ) * q₁(κ, μ), sqrt(κ) * q₂(κ, μ))
end

function f₂(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) + μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) - Γ(κ, μ) * term
end

function Intθ_2(κ, μ)
    return (sqrt(κ) * q₂(κ, μ), pi - sqrt(κ) * q₁(κ, μ))
end

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

function f₄(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) - μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) + Γ(κ, μ) * term
end

function Intθ_4(κ, μ)
    return (pi - sqrt(κ) * q₁(κ, μ), pi + sqrt(κ) * q₂(κ, μ))
end

### Derivatives
function f₁θ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    Γv = Γ(κ, μ)

    d_main = cos(twice(θ))

    sa = sin(a)
    ca = cos(a)
    sb = sin(b)
    cb = cos(b)

    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)
    d_term = (-twice(ca) * inv_sa3) + (twice(μ) * cb * inv_sb3)

    return d_main - Γv * d_term
end





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


# # compare analytic vs automatic derivative
# κ = 0.5
# μ = 0.02
# θ = 0.7
# ad = ForwardDiff.derivative(θ -> f₁(θ, κ, μ), θ)
# an = f₁θ(θ, κ, μ)

# println("analytic  f₁θ = ", an)
# println("ForwardDiff = ", ad)
# println("difference  = ", ad - an)

# bounds for the box
μT = 0.0886765
κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, μT

# make an interior grid (with margin)
interior_linspace(lo, hi, n; margin=0.02) = begin
    a = lo + margin * (hi - lo)
    b = hi - margin * (hi - lo)
    range(a, b; length=n)
end

# Example: make all grid pairs as a vector of tuples (κ, μ)
nk, nm = 20, 20
κs = collect(interior_linspace(κ_lo, κ_hi, nk))  # κ-grid
μs = collect(interior_linspace(μ_lo, μ_hi, nm))  # μ-grid

grid = [(κ, μ) for κ in κs, μ in μs]  # nk×nm matrix of pairs
roots_grid = zeros(nk * nm)  # storage for roots at each (κ, μ)
pairs = vec(grid)  # length = nk*nm

function brackets_signchanges(f, a::Real, b::Real; n::Integer=400)
    xs = collect(range(a, b; length=n + 1))
    br = Tuple{Float64,Float64}[]

    # evaluate first finite value
    i = 1
    fa = f(xs[i])
    while i <= n && !isfinite(fa)
        i += 1
        fa = f(xs[i])
    end
    i > n && return br  # nothing finite

    sgn(x) = x > 0 ? 1 : x < 0 ? -1 : 0
    sa = sgn(fa)

    for j in i:n
        x1, x2 = xs[j], xs[j+1]
        fx = f(x2)
        if !isfinite(fx)
            continue  # skip nonfinite; keep previous finite fa/sa
        end
        sb = sgn(fx)

        if sb == 0
            push!(br, (x2, x2))              # exact root at node
        elseif sa == 0
            push!(br, (x1, x1))              # exact root at previous node
        elseif sa != sb
            push!(br, (x1, x2))              # sign flip bracket
        end

        fa, sa = fx, sb                      # advance
    end
    return br
end

# Brackets for f₁ at fixed (κ, μ) over θ ∈ (-√κ q₁, √κ q₂), trimmed by `margin`
function brackets_f₁(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = -sqrt(κ) * q₁(κ, μ)
    θmax = sqrt(κ) * q₂(κ, μ)
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₁(θ, κ, μ), a, b; n=n)
end


function brackets_f₂(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_2(κ, μ)[1]
    θmax = Intθ_2(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₂(θ, κ, μ), a, b; n=n)
end

function brackets_f₃(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_3(κ, μ)[1]
    θmax = Intθ_3(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₃(θ, κ, μ), a, b; n=n)
end

function brackets_f₄(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_4(κ, μ)[1]
    θmax = Intθ_4(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₄(θ, κ, μ), a, b; n=n)
end

function pairs_to_roots(pairs; n=400)
    roots = zeros(length(pairs))
    for i in 1:length(pairs)
        κ, μ = pairs[i]
        br = brackets_f₁(κ, μ; n=n)
        println("i=$i, (κ,μ)=($κ,$μ), brackets = ", br)
        if length(br) == 0
            roots[i] = NaN  # no root found
        else
            rts = find_zero(θ -> f₁(θ, κ, μ), br[1])
            roots[i] = rts[1]  # take first root
        end
    end
    return roots
end

roots = pairs_to_roots(pairs; n=400)



### Plotting f1
# xs = first.(pairs)           # κ
# ys = last.(pairs)            # μ

# # collect sorted unique axes
# κs = sort(unique(first.(pairs)))
# μs = sort(unique(last.(pairs)))
# nk, nm = length(κs), length(μs)

# # fill θ*(κ_i, μ_j) matrix
# θgrid = fill(NaN, nk, nm)
# for (i, (κ, μ)) in enumerate(pairs)
#     ii = searchsortedfirst(κs, κ)
#     jj = searchsortedfirst(μs, μ)
#     θgrid[ii, jj] = roots[i]
# end

# Plots.surface(μs, κs, θgrid;
#     xlabel="μ", ylabel="κ", zlabel="θ*",
#     title="Root surface θ*(κ, μ)",
#     legend=false,
#     color=:cork, cbar_title="θ*",
#     linewidth=0,
#     fillalpha=0.95,
#     camera=(30, 60),
#     grid=:on, framestyle=:box)



# ### A single root Example 
κ, μ = 0.5, 0.02
# ### Estimate by sign changes
# for i in range(Intθ(κ, μ)[1],Intθ(κ, μ)[2]; length=400)
#     println("f₁($i) = ", f₁(i, κ, μ))
# end


# br = brackets_f₁(κ, μ; n=800)
# f(θ) = f₁(θ, κ, μ)
# rts = refine_brackets(f, br)
# println(rts)
# find_zero(θ -> f₁(θ, κ, μ), br[1])


br = brackets_f₂(κ, μ; n=800)

# rts = refine_brackets(f, br)
# println(rts)
# find_zero(θ -> f₁(θ, κ, μ), br[1])
