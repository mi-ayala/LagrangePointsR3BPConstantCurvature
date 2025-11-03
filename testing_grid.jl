# Plain numeric rectangle for (κ, μ)
struct Rect
    κlo::Float64
    κhi::Float64
    μlo::Float64
    μhi::Float64
    ok::Bool      # default: false until validated
    level::Int    # refinement level (0 for seeds)
end

# convenience constructors
Rect(κlo, κhi, μlo, μhi) = Rect(κlo, κhi, μlo, μhi, false, 0)

# Regular nx×ny grid inside a parent rectangle → vector of Rects (level can be set)
function grid_rect(parent::Rect; nx::Int, ny::Int, level::Int=1)
    @assert nx > 0 && ny > 0
    κstep = (parent.κhi - parent.κlo) / nx
    μstep = (parent.μhi - parent.μlo) / ny

    cells = Vector{Rect}(undef, nx * ny)
    t = 0
    for i in 0:nx-1, j in 0:ny-1
        t += 1
        cells[t] = Rect(parent.κlo + i * κstep, parent.κlo + (i + 1) * κstep,
            parent.μlo + j * μstep, parent.μlo + (j + 1) * μstep,
            false, level)
    end
    return cells
end

μT = 0.0886765
κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, μT
R = Rect(κ_lo, κ_hi, μ_lo, μ_hi)       # whole domain, level 0


cells = grid_rect(R; nx=4, ny=3, level=1)

println(length(cells), " cells")
println(cells[1])        # first cell
println(cells[end])      # last cell
