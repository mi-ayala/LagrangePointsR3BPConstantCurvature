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

### Functions and interval
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
    y = tan(r * q₂(κ, μ))
    x = one(μ) + μ
    return atan(y, x)
end


# left endpoint of the above interval
function θ4_min(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    πT = π * one(κ)                    # π in the promoted type

    # sqrt( 4Γ / (2κ q2^2 + π^2 Γ μ) )
    oneG = one(G)
    two = oneG + oneG
    four = two + two
    denom = two * κ * q2^2 + (πT^2) * G * μ
    rad = (four * G) / denom
    s = sqrt(rad)

    ones = one(s)                      # 1 in the same type as s
    return πT - (ones - s) * r * q2
end

# right endpoint of the interval
function θ4_max(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    πT = π * one(κ)

    # sqrt( 4Γ μ / (2κ q1^2 + π^2 Γ) )
    oneG = one(G)
    two = oneG + oneG
    four = two + two
    denom = two * κ * q1^2 + (πT^2) * G
    rad = (four * G * μ) / denom
    s = sqrt(rad)

    ones = one(s)
    return πT + (ones - (q1 / q2) * s) * r * q1
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


# f₄θ(θ, κ, μ) = ∂f₄/∂θ
function f₄θ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    Γv = Γ(κ, μ)

    a = θ + r * q1
    b = θ - r * q2

    d_main = cos(twice(θ))

    sa = sin(a)
    ca = cos(a)
    sb = sin(b)
    cb = cos(b)

    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)

    d_term = (-twice(ca) * inv_sa3) + (twice(μ) * cb * inv_sb3)

    return d_main + Γv * d_term   # note the + (for f₁θ it was -)
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

    μT = 0.0886765
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

"""
    plot_validated_boxes(boxes;
        mode=:validated,                        # :validated | :invalid | :both
        color_valid=:green, color_invalid=:orange,
        alpha=0.4, edgecolor=:black, lw=0.5,
        markersize=4, title="Boxes",
        show_bounding_box=true)

Plot boxes of the form (ok::Bool, κ, μ, rκ, rμ).
`mode=:validated` (default) plots only ok=true,
`mode=:invalid` plots only ok=false,
`mode=:both` plots both (different colors).

Adds the bounding region κ ∈ [0, (0.5π)^2], μ ∈ [0, μT≈0.0886765].
"""
function plot_validated_boxes(boxes;
    mode::Symbol=:validated,
    color_valid=:green, color_invalid=:orange,
    alpha::Real=0.4, edgecolor=:black, lw::Real=3,
    markersize::Real=4, title::AbstractString="Boxes",
    show_bounding_box::Bool=true)

    μT = 0.0886765
    κ_lo, κ_hi = 0.0, (0.5 * π)^2
    μ_lo, μ_hi = 0.0, μT

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        title=title, legend=false,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    # helper to draw one box
    draw_box!(κ, μ, rκ, rμ; col) = begin
        xs = [κ - rκ, κ + rκ, κ + rκ, κ - rκ, κ - rκ]
        ys = [μ - rμ, μ - rμ, μ + rμ, μ + rμ, μ - rμ]
        Plots.plot!(plt, xs, ys; color=col, alpha=alpha, linecolor=edgecolor, lw=lw)
        # Plots.scatter!(plt, [κ], [μ]; color=edgecolor, ms=markersize)
    end

    if mode == :validated || mode == :both
        for (ok, κ, μ, rκ, rμ) in boxes
            ok || continue
            draw_box!(κ, μ, rκ, rμ; col=color_valid)
        end
    end
    if mode == :invalid || mode == :both
        for (ok, κ, μ, rκ, rμ) in boxes
            ok && continue
            draw_box!(κ, μ, rκ, rμ; col=color_invalid)
        end
    end

    if show_bounding_box
        Xb = [κ_lo, κ_hi, κ_hi, κ_lo, κ_lo]
        Yb = [μ_lo, μ_lo, μ_hi, μ_hi, μ_lo]
        Plots.plot!(plt, Xb, Yb; color=:transparent, linecolor=:red, lw=3)
    end

    display(plt)
    return plt
end

"""
    plot_validated_boxes(
        boxes;
        mode=:validated,                    # :validated | :invalid | :both
        color_valid=:green, color_invalid=:red,
        alpha=0.4, edgecolor=:black, lw=0.5,
        markersize=3, title="Boxes",
        show_bounding_box=true,
        # NEW:
        y_exaggerate::Real=1.0,             # multiply rμ for drawing only
        min_frac_h::Real=0.0,               # enforce min visual height = min_frac_h*(μ_hi-μ_lo)
        plot_centers_invalid::Bool=true,    # scatter centers of invalid boxes
        plot_centers_valid::Bool=false,     # (optional) centers of valid boxes
        zoom::Union{Nothing,Tuple{Real,Real,Real,Real}}=nothing # (κlo, κhi, μlo, μhi)
    )

Visualize boxes (ok, κ, μ, rκ, rμ). To make tiny boxes visible, set
`y_exaggerate > 1` or `min_frac_h > 0`. These affect drawing only.
"""

"""
    plot_validated_boxes(
        boxes;
        mode=:validated,                    # :validated | :invalid | :both
        color_valid=:green, color_invalid=:red,
        alpha=0.4, edgecolor=:black, lw=0.5,
        markersize=3, title="Boxes",
        show_bounding_box=true,
        # NEW:
        y_exaggerate::Real=1.0,             # multiply rμ for drawing only
        min_frac_h::Real=0.0,               # enforce min visual height = min_frac_h*(μ_hi-μ_lo)
        plot_centers_invalid::Bool=true,    # scatter centers of invalid boxes
        plot_centers_valid::Bool=false,     # (optional) centers of valid boxes
        zoom::Union{Nothing,Tuple{Real,Real,Real,Real}}=nothing # (κlo, κhi, μlo, μhi)
    )

Visualize boxes (ok, κ, μ, rκ, rμ). To make tiny boxes visible, set
`y_exaggerate > 1` or `min_frac_h > 0`. These affect drawing only.
"""
function plot_validated_boxes(boxes;
    mode::Symbol=:validated,
    color_valid=:green, color_invalid=:red,
    alpha::Real=0.4, edgecolor=:black, lw::Real=0.5,
    markersize::Real=1, title::AbstractString="Boxes",
    show_bounding_box::Bool=true,
    y_exaggerate::Real=1.0,
    min_frac_h::Real=0.0,
    plot_centers_invalid::Bool=true,
    plot_centers_valid::Bool=false,
    zoom::Union{Nothing,Tuple{Real,Real,Real,Real}}=nothing)

    μT = 0.0886765
    κ_lo_full, κ_hi_full = 0.0, (0.5 * π)^2
    μ_lo_full, μ_hi_full = 0.0, μT

    # axes limits
    κ_lo, κ_hi, μ_lo, μ_hi = isnothing(zoom) ? (κ_lo_full, κ_hi_full, μ_lo_full, μ_hi_full) : zoom
    μ_span = μ_hi - μ_lo
    min_h = max(0.0, min_frac_h) * μ_span

    plt = Plots.plot(; xlabel="κ", ylabel="μ", title=title, legend=false,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    # helper: draw one rectangle with visual tweaks on height only
    function draw_box!(κ, μ, rκ, rμ; col)
        rμ_vis = max(rμ * y_exaggerate, min_h / 2)   # visual-only height
        xs = [κ - rκ, κ + rκ, κ + rκ, κ - rκ, κ - rκ]
        ys = [μ - rμ_vis, μ - rμ_vis, μ + rμ_vis, μ + rμ_vis, μ - rμ_vis]
        Plots.plot!(plt, xs, ys; color=col, alpha=alpha, linecolor=edgecolor, lw=lw)
    end

    # optionally collect centers for scatter overlays
    κ_invalid = Float64[]
    μ_invalid = Float64[]
    κ_valid = Float64[]
    μ_valid = Float64[]

    # draw according to mode
    if mode == :validated || mode == :both
        for (ok, κ, μ, rκ, rμ) in boxes
            ok || continue
            draw_box!(κ, μ, rκ, rμ; col=color_valid)
            plot_centers_valid && (push!(κ_valid, κ);
            push!(μ_valid, μ))
        end
    end
    if mode == :invalid || mode == :both
        for (ok, κ, μ, rκ, rμ) in boxes
            ok && continue
            draw_box!(κ, μ, rκ, rμ; col=color_invalid)
            plot_centers_invalid && (push!(κ_invalid, κ);
            push!(μ_invalid, μ))
        end
    end

    # overlay centers (helps when rectangles are too thin)
    if plot_centers_invalid && !isempty(κ_invalid)
        Plots.scatter!(plt, κ_invalid, μ_invalid; color=color_invalid, ms=markersize)
    end
    if plot_centers_valid && !isempty(κ_valid)
        Plots.scatter!(plt, κ_valid, μ_valid; color=color_valid, ms=markersize)
    end

    # outer bounding box (full region), unless zooming makes it confusing
    if show_bounding_box
        Xb = [κ_lo_full, κ_hi_full, κ_hi_full, κ_lo_full, κ_lo_full]
        Yb = [μ_lo_full, μ_lo_full, μ_hi_full, μ_hi_full, μ_lo_full]
        Plots.plot!(plt, Xb, Yb; color=:transparent, linecolor=:red, lw=3)
    end

    display(plt)
    return plt
end


# helpers (safe for Reals and Intervals)
@inline twice(x) = x + x
@inline csc(x) = one(x) / sin(x)
@inline powabs3(x) = abs(x)^3



# ---- Notable functions ----
# λ₁(θ) = -[ cos(2θ) + 2Γ(  cos(θ+√κ q₁)/|sin(θ+√κ q₁)|^3  + μ cos(θ-√κ q₂)/|sin(θ-√κ q₂)|^3 ) ]
function λ₁(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    s1 = sin(θ + r * q1)
    s2 = sin(θ - r * q2)
    t1 = cos(θ + r * q1) / powabs3(s1)
    t2 = cos(θ - r * q2) / powabs3(s2)

    return -(cos(twice(θ)) + twice(G) * (t1 + μ * t2))
end

# λ₂(θ) = -Γ sinθ [ sin(√κ)/|sin(θ+√κ q₁)|^3  - μ sin(√κ)/|sin(θ-√κ q₂)|^3 ]
function λ₂(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    s1 = sin(θ + r * q1)
    s2 = sin(θ - r * q2)
    sκ = sin(r)

    return -G * sin(θ) * (sκ / powabs3(s1) - μ * (sκ / powabs3(s2)))
end

# det(θ) = λ₁(θ) * λ₂(θ)
detλ(θ, κ, μ) = λ₁(θ, κ, μ) * λ₂(θ, κ, μ)

# c₂(θ) = 2 + cos(2θ) - Γ csc(θ) [ sin(2θ+√κ q₁)/|sin(θ+√κ q₁)|^3 - μ sin(2θ-√κ q₂)/|sin(θ-√κ q₂)|^3 ]
function c₂(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    s1 = sin(θ + r * q1)
    s2 = sin(θ - r * q2)
    a = sin(twice(θ) + r * q1) / powabs3(s1)
    b = sin(twice(θ) - r * q2) / powabs3(s2)

    two = one(θ) + one(θ)      # “2” with the same type as θ (Float64, Interval, BigFloat,…)
    return two + cos(twice(θ)) - G * csc(θ) * (a - μ * b)
end

# Δ(θ) = c₂(θ)^2 - 4*det(θ)/sin^2(θ)
function Δ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    c2 = c₂(θ, κ, μ)
    four = 4 * one(c2)
    sθ = sin(θ)
    two = one(θ) + one(θ)
    four = two + two
    return c2^2 - (four * detλ(θ, κ, μ)) / (sθ^2)
end


