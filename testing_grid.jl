# --- minimal Rect + grid helper (as used before) ---
struct Rect
    κlo::Float64; κhi::Float64
    μlo::Float64; μhi::Float64
    ok::Bool
    level::Int
end
Rect(κlo, κhi, μlo, μhi) = Rect(κlo, κhi, μlo, μhi, false, 0)

function grid_rect(parent::Rect; nx::Int, ny::Int, level::Int=1)
    @assert nx > 0 && ny > 0
    κstep = (parent.κhi - parent.κlo) / nx
    μstep = (parent.μhi - parent.μlo) / ny
    cells = Vector{Rect}(undef, nx*ny)
    t = 0
    for i in 0:nx-1, j in 0:ny-1
        t += 1
        cells[t] = Rect(parent.κlo + i*κstep, parent.κlo + (i+1)*κstep,
                        parent.μlo + j*μstep, parent.μlo + (j+1)*μstep,
                        false, level)
    end
    return cells
end

# --- build inner rectangle with fixed margins, then a 3×3 grid inside it ---
"""
    inner_3x3(κ_lo, κ_hi, μ_lo, μ_hi; margin=0.002)

Return `(inner::Rect, cells::Vector{Rect})` where `inner` is the domain
shrunk by `margin` on all sides, and `cells` is a regular 3×3 grid (level=1)
covering `inner`.
"""
function inner_3x3(κ_lo, κ_hi, μ_lo, μ_hi; margin=0.002)
    # clamp margins if user passes something too large
    Lκ = κ_hi - κ_lo
    Lμ = μ_hi - μ_lo
    mκ = min(margin, Lκ/2 - eps())   # ensure positive inner width
    mμ = min(margin, Lμ/2 - eps())   # ensure positive inner height

    inner = Rect(κ_lo + mκ, κ_hi - mκ, μ_lo + mμ, μ_hi - mμ, false, 0)
    cells = grid_rect(inner; nx=3, ny=3, level=1)
    return inner, cells
end

# --- quick smoke test ---
# Domain defaults: κ ∈ [0, (π/2)^2], μ ∈ [0, μT]
κ_lo, κ_hi = 0.0, (0.5*pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

inner, cells = inner_3x3(κ_lo, κ_hi, μ_lo, μ_hi; margin=0.002)
println("Inner rect: ", inner)
println("3×3 cells = ", length(cells))
println("First cell: ", cells[1])
println("Last  cell: ", cells[end])


cells[8]