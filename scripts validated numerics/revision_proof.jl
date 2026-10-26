using RC3BP
using RadiiPolynomial
using Roots


# FOR MIGUEL: Perform a CAP to show that the function $g_{\kappa,\mu}$ in Eq. (6.5) (already coded in Github) has exactly 4 zeros in the interval
# $[a_{\kappa,\mu},b_{\kappa,\mu}]$ for values of $\kappa, \mu$ close to $\kappa=0.575045$ and $\mu=0.20591$ (provide bounds). The domain of $g_{\kappa,\mu}$ is approximately  $[a_{\kappa,\mu},b_{\kappa,\mu}]=[0.394367,1.85444]$ and
# the roots are approximately $x_1=0.872479$, $x_2=1.73155$, $x_3=1.79198$, $x_4=1.85231$. It would also be useful, although less important, 
# to calculate the number of collinear RE for this values of $\kappa, \mu$ and to run the stability routines. This in order to count the `maximum number of RE' and
# also to check consistency with Theorem 4.1.


κ = 0.575045
μ = 0.20591

Idomain = [0.394367, 1.85444]

x1 = find_zero(θ -> g(θ, κ, μ), 0.872479; atol=1e-12)
x2 = find_zero(θ -> g(θ, κ, μ), 1.73155; atol=1e-12)
x3 = find_zero(θ -> g(θ, κ, μ), 1.79198; atol=1e-12)
x4 = find_zero(θ -> g(θ, κ, μ), 1.85231; atol=1e-12)

Iκ = interval(κ, 1e-6; format=:midpoint)
Iμ = interval(μ, 1e-6; format=:midpoint)

#### Plot numerics
using CairoMakie
fig = Figure()
ax = Axis(
    fig[1, 1];
    xlabel="θ",
    ylabel="g(θ, κ, μ)",
    title="Function g(θ, κ, μ) for κ=$κ and μ=$μ",
    xgridvisible=true,
    ygridvisible=true,
)
θs = range(1.7, Idomain[2]; length=1000)
ys = [g(θ, κ, μ) for θ in θs]
lines!(ax, θs, ys; label="g(θ, κ, μ)")
fig


### Validate zeros
x1_int = validate_root_sign_mono(
    g,
    gθ,
    x1,
    interval(x1, 1e-1; format=:midpoint),
    Iκ,
    Iμ
)

x2_int = validate_root_sign_mono(
    g,
    gθ,
    x2,
    interval(x2, 0.5; format=:midpoint),
    Iκ,
    Iμ
)

x3_int = validate_root_sign_mono(
    g,
    gθ,
    x3,
    interval(x3, 1e-2; format=:midpoint),
    Iκ,
    Iμ
)

x4_int = validate_root_sign_mono(
    g,
    gθ,
    x4,
    interval(x4, 1e-3; format=:midpoint),
    Iκ,
    Iμ
)

println("x₁ validated: ", x1_int[1], ", root interval: ", x1_int[2].Iexist)
println("x₂ validated: ", x2_int[1], ", root interval: ", x2_int[2].Iexist)
println("x₃ validated: ", x3_int[1], ", root interval: ", x3_int[2].Iexist)
println("x₄ validated: ", x4_int[1], ", root interval: ", x4_int[2].Iexist)


### Check that g is nonzero in the complement of the intervals around the roots

Icomp1 = interval(Idomain[1], inf(x1_int[2].J))
Icomp2 = interval(sup(x1_int[2].J), inf(x2_int[2].J))
Icomp3 = interval(sup(x2_int[2].J), inf(x3_int[2].J))
Icomp4 = interval(sup(x3_int[2].J), inf(x4_int[2].J))
Icomp5 = interval(sup(x4_int[2].J), Idomain[2])

Ig1 = eval_subdiv(g, Icomp1, Iκ, Iμ; ntiles=2)
Ig2 = eval_subdiv(g, Icomp2, Iκ, Iμ; ntiles=5000)
Ig3 = eval_subdiv(g, Icomp3, Iκ, Iμ; ntiles=500)
Ig4 = eval_subdiv(g, Icomp4, Iκ, Iμ; ntiles=5000)
Ig5 = eval_subdiv(g, Icomp5, Iκ, Iμ; ntiles=500)

println("g on ", Icomp1, " = ", Ig1, "; excludes zero: ", nonzero(Ig1))
println("g on ", Icomp2, " = ", Ig2, "; excludes zero: ", nonzero(Ig2))
println("g on ", Icomp3, " = ", Ig3, "; excludes zero: ", nonzero(Ig3))
println("g on ", Icomp4, " = ", Ig4, "; excludes zero: ", nonzero(Ig4))
println("g on ", Icomp5, " = ", Ig5, "; excludes zero: ", nonzero(Ig5))
