function brackets_signchanges(f, a::Real, b::Real; n::Integer=200)
    xs = collect(range(a, b; length=n + 1))
    br = Tuple{Float64,Float64}[]

    i = 1
    fa = f(xs[i])
    while i <= n && !isfinite(fa)
        i += 1
        fa = f(xs[i])
    end
    i > n && return br

    sgn(x) = x > 0 ? 1 : x < 0 ? -1 : 0
    sa = sgn(fa)

    for j in i:n
        x1, x2 = xs[j], xs[j+1]
        fx = f(x2)
        if !isfinite(fx)
            continue
        end
        sb = sgn(fx)

        if sb == 0
            push!(br, (x2, x2))              
        elseif sa == 0
            push!(br, (x1, x1))              
        elseif sa != sb
            push!(br, (x1, x2))             
        end

        fa, sa = fx, sb                     
    end
    return br
end

function classify_equilibrium_f(θ, κ, μ; verbose::Bool=false)
    det_val = detλ(θ, κ, μ)
    lam2_val = λ₂(θ, κ, μ)
    Δ_val = Δ(θ, κ, μ)

    label =
        det_val < 0 ? :center_saddle_unstable :
        lam2_val > 0 ? :lyapunov_stable :
        Δ_val > 0 ? :elliptic :
        Δ_val < 0 ? :complex_saddle_unstable :
        :unknown

    if verbose
        println("θ=$θ κ=$κ μ=$μ → label=$label, det=$det_val, λ₂=$lam2_val, Δ=$Δ_val")
    end

    return label
end

function classify_center_equilibrium_num(x0::Real, κ::Real, μ::Real; verbose::Bool=false)
    x0 = float(x0)
    κ = float(κ)
    μ = float(μ)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)

    Λ = cbrt(sin(r * q1) / (μ * sin(r * q2)))

    arg = Λ * sin(x0)
    if !(isfinite(arg) && -1.0 ≤ arg ≤ 1.0)
        verbose && println("asin argument not in [-1,1] or not finite: ", arg)
        z = NaN
        return false, :unknown, z, z, z, z, x0, Λ
    end

    d1 = asin(arg)
    d2 = x0

    sr = sin(r)
    sr4 = sr^4

    s_rq1 = sin(r * q1)
    s_rq2 = sin(r * q2)

    pref = (3.0 * s_rq1 * s_rq2) / sr4

    B1 = s_rq2 * cos(d1) + s_rq1 * cos(d2)

    B2 = (1.0 + 2.0 * cos(2.0 * d2)) * s_rq1 * cos(d1) +
         (1.0 + 2.0 * cos(2.0 * d1)) * s_rq2 * cos(d2)

    det = pref * B1 * B2

    tr = (1.0 / (sr^2)) * (
        (1.0 - 4.0 * cos(d2)^2) * s_rq1^2
        -
        6.0 * cos(d1) * cos(d2) * s_rq1 * s_rq2
        +
        (1.0 - 4.0 * cos(d1)^2) * s_rq2^2
    )

    m = -(sinsh = (sin(d1)^2 * sin(d2)^2)) /
        ((cos(r) - cos(d1 + d2)) * (cos(r) - cos(d1 - d2)))

    C = m - 4.0 * det

    verbose && println("det=$det  tr=$tr  C=$C")

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

