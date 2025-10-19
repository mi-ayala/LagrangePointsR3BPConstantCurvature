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

function θ1_max(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    y = tan(r * q₂(κ, μ))     # numerator
    x = one(μ) + μ            # denominator (kept as same type)
    return atan(y, x)
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


### Plotting function 


"""
    plot_validated_boxes(boxes;
        color=:green, alpha=0.4, edgecolor=:black,
        markersize=4, title="Validated boxes",
        show_bounding_box=true)

Plot validated boxes of the form (ok::Bool, κ, μ, rκ, rμ),
centered at (κ, μ) with half-width rκ and half-height rμ.

Adds the bounding parameter region:
    κ ∈ [0, (0.5π)^2],  μ ∈ [0, μT≈0.0886765].
"""
function plot_validated_boxes(boxes;
        color=:green, alpha=0.4, edgecolor=:black,
        markersize=4, title="Validated boxes",
        show_bounding_box=true)

    μT   = 0.0886765
    κ_lo, κ_hi = 0.0, (0.5 * π)^2
    μ_lo, μ_hi = 0.0, μT

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
     title=title, legend=false,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    # draw all validated boxes
    for (ok, κ, μ, rκ, rμ) in boxes
        ok || continue
        xs = [κ - rκ, κ + rκ, κ + rκ, κ - rκ, κ - rκ]
        ys = [μ - rμ, μ - rμ, μ + rμ, μ + rμ, μ - rμ]
        Plots.plot!(plt, xs, ys; color=color, alpha=alpha, linecolor=edgecolor, lw=0.5)
        # Plots.scatter!(plt, [κ], [μ]; color=edgecolor, ms=markersize)
    end

    # bounding box
    if show_bounding_box
        Xb = [κ_lo, κ_hi, κ_hi, κ_lo, κ_lo]
        Yb = [μ_lo, μ_lo, μ_hi, μ_hi, μ_lo]
        Plots.plot!(plt, Xb, Yb; color=:transparent, linecolor=:red, lw=3, label="Bounding region")
    end

    Plots.display(plt)
    return plt
end

