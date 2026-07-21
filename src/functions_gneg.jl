function q₁neg(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(-κ)
    y = μ * sinh(twice(r))
    x = one(μ) + μ * cosh(twice(r))
    return (one(r) / twice(r)) * atanh(y / x)
end

function q₂neg(κ, μ)
    κ, μ = promote(κ, μ)
    return one(κ) - q₁neg(κ, μ)
end

function Γneg(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(-κ)
    num = sqr(sinh(r)) * sinh(r) * cosh(r)
    den = sqrt(μ * μ + twice(μ * cosh(twice(r))) + one(μ))
    return num / den
end

function Λneg(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(-κ)
    q1 = q₁neg(κ, μ)
    q2 = q₂neg(κ, μ)

    s1 = sinh((-r) * q1)
    s2 = sinh((-r) * q2)

    return cbrt_t(s1 / (μ * s2))
end

function gneg(x, κ, μ)
    x, κ, μ = promote(x, κ, μ)

    r = sqrt(-κ)
    q1 = q₁neg(κ, μ)
    q2 = q₂neg(κ, μ)
    G = Γneg(κ, μ)
    L = Λneg(κ, μ)

    s1 = sinh((-r) * q1)
    s2 = sinh((-r) * q2)

    shx = sinh(x)
    chx = cosh(x)

    inner = one(shx) + (L * L) * (shx * shx)
    rt = sqrt(inner)

    bracket = s1 * chx + s2 * rt

    pref = (L * L * L) * s2 * (shx * shx * shx)
    tail = G * (sinh(-r) * sinh(-r))

    return pref * bracket - tail
end

function classify_center_equilibrium_num_neg(x0::Real, κ::Real, μ::Real; verbose::Bool=false)
    x0 = float(x0)
    κ = float(κ)
    μ = float(μ)

    r = sqrt(-κ)
    q1 = q₁neg(κ, μ)
    q2 = q₂neg(κ, μ)

    Λ = cbrt(cosh((-r) * q2) / cosh((-r) * q1))

    d1 = asinh(Λ * sinh(x0))
    d2 = x0

    sh_r = sinh(-r)
    sh_r2 = sh_r * sh_r
    sh_r4 = sh_r2 * sh_r2

    sh_rq1 = sinh((-r) * q1)
    sh_rq2 = sinh((-r) * q2)

    ch_r = cosh(-r)

    ch_d1 = cosh(d1)
    ch_d2 = cosh(d2)

    pref = (3.0 * sh_rq1 * sh_rq2) / sh_r4

    B1 = sh_rq2 * ch_d1 + sh_rq1 * ch_d2
    B2 = (1.0 + 2.0 * cosh(2.0 * d2)) * sh_rq1 * ch_d1 +
         (1.0 + 2.0 * cosh(2.0 * d1)) * sh_rq2 * ch_d2

    det = pref * B1 * B2


    tr = (1.0 / sh_r2) * (
        (1.0 - 4.0 * (ch_d2 * ch_d2)) * (sh_rq1 * sh_rq1)
        -
        6.0 * ch_d1 * ch_d2 * sh_rq1 * sh_rq2
        +
        (1.0 - 4.0 * (ch_d1 * ch_d1)) * (sh_rq2 * sh_rq2)
    )

    sh_d1 = sinh(d1)
    sh_d2 = sinh(d2)
    num_m = (sh_d1 * sh_d1) * (sh_d2 * sh_d2)
    den_m = (ch_r - cosh(d1 + d2)) * (ch_r - cosh(d1 - d2))
    m = -(num_m / den_m)

    C = m - 4.0 * det

    verbose && println("Λ=$Λ  d1=$d1  d2=$d2  det=$det  tr=$tr  m=$m  C=$C")

    if !(isfinite(det) && isfinite(tr) && isfinite(C))
        return :unknown
    end

    if det < 0
        return :center_saddle_unstable
    elseif tr > 0
        return :lyapunov_stable
    elseif C > 0
        return :elliptic
    else
        return :complex_saddle_unstable
    end
end

function Intθ_gneg(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(-κ)
    L = Λneg(κ, μ)

    num = sinh(-r)
    den = L + cosh(-r)

    return (atanh(num / den), 10.0)
end
