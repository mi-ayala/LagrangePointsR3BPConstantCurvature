using RC3BP
using RadiiPolynomial
using Roots


κ = 0.575045
μ = 0.20591

Idomain = [0.394367, 1.85444]

Iκ = interval(κ, 1e-5; format=:midpoint)
Iμ = interval(μ, 1e-5; format=:midpoint)

x1 = find_zero(θ -> g(θ, κ, μ), 0.872479; atol=1e-12)
x2 = find_zero(θ -> g(θ, κ, μ), 1.73155; atol=1e-12)
x3 = find_zero(θ -> g(θ, κ, μ), 1.79198; atol=1e-12)
x4 = find_zero(θ -> g(θ, κ, μ), 1.85231; atol=1e-12)

I1 = [Idomain[1], x1]
I2 = [x1, x2]
I3 = [x2, x3]
I4 = [x3, x4]
I5 = [x4, Idomain[2]]


xmin_I1 = I1[1]
xmax_I2 = find_zero(θ -> gθ(θ, κ, μ), (I2[1], I2[2]))
xmin_I3 = find_zero(θ -> gθ(θ, κ, μ), (I3[1], I3[2]))
xmax_I4 = find_zero(θ -> gθ(θ, κ, μ), (I4[1], I4[2]))
xmin_I5 = I5[2]

min_g_I1 = g(xmin_I1, κ, μ)
max_g_I2 = g(xmax_I2, κ, μ)
min_g_I3 = g(xmin_I3, κ, μ)
max_g_I4 = g(xmax_I4, κ, μ)
min_g_I5 = g(xmin_I5, κ, μ)

println("min g on I1: $min_g_I1 at θ = $xmin_I1; sign = $(sign(min_g_I1))")
println("max g on I2: $max_g_I2 at θ = $xmax_I2; sign = $(sign(max_g_I2))")
println("min g on I3: $min_g_I3 at θ = $xmin_I3; sign = $(sign(min_g_I3))")
println("max g on I4: $max_g_I4 at θ = $xmax_I4; sign = $(sign(max_g_I4))")
println("min g on I5: $min_g_I5 at θ = $xmin_I5; sign = $(sign(min_g_I5))")


Ig1 = eval_subdiv(g, interval(xmin_I1), Iκ, Iμ; ntiles=50)
Ig2 = eval_subdiv(g, interval(xmax_I2), Iκ, Iμ; ntiles=50)
Ig3 = eval_subdiv(g, interval(xmin_I3), Iκ, Iμ; ntiles=50)
Ig4 = eval_subdiv(g, interval(xmax_I4), Iκ, Iμ; ntiles=50)
Ig5 = eval_subdiv(g, interval(xmin_I5), Iκ, Iμ; ntiles=500)

println("interval evaluation at min on I1: $Ig1; sign = $(sign_tag(Ig1))")
println("interval evaluation at max on I2: $Ig2; sign = $(sign_tag(Ig2))")
println("interval evaluation at min on I3: $Ig3; sign = $(sign_tag(Ig3))")
println("interval evaluation at max on I4: $Ig4; sign = $(sign_tag(Ig4))")
println("interval evaluation at min on I5: $Ig5; sign = $(sign_tag(Ig5))")
