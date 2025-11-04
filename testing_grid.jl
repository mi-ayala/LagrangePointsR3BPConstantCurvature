using IntervalArithmetic
using Plots
using PlotlyJS
using Roots


struct Rect
    κlo::Float64
    κhi::Float64
    μlo::Float64
    μhi::Float64
    ok::Bool
    level::Int
end
Rect(κlo, κhi, μlo, μhi) = Rect(κlo, κhi, μlo, μhi, false, 0)

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
    mκ = min(margin, Lκ / 2 - eps())   # ensure positive inner width
    mμ = min(margin, Lμ / 2 - eps())   # ensure positive inner height

    inner = Rect(κ_lo + mκ, κ_hi - mκ, μ_lo + mμ, μ_hi - mμ, false, 0)
    cells = grid_rect(inner; nx=3, ny=3, level=1)
    return inner, cells
end

# --- quick smoke test ---
# Domain defaults: κ ∈ [0, (π/2)^2], μ ∈ [0, μT]
κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

inner, cells = inner_3x3(κ_lo, κ_hi, μ_lo, μ_hi; margin=0.002)
println("Inner rect: ", inner)
println("3×3 cells = ", length(cells))
println("First cell: ", cells[1])
println("Last  cell: ", cells[end])


cells[5]


include("functions.jl")
include("validation_functions.jl")
"""
    cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi;
        nx::Int, ny::Int, bracket_n::Int=400, verbose::Bool=true)

Partition the given rectangle into an `nx × ny` regular grid (no margins).
For each cell, take (rκ, rμ) = (Δκ/2, Δμ/2), bracket & refine a root at the
cell center, then validate with `validate_f1_root`.

Returns:
  boxes :: Vector{Tuple{Bool,Float64,Float64,Float64,Float64}}
           (ok, κc, μc, rκ, rμ) for validated cells
  stats :: NamedTuple with {Δκ, Δμ, rκ, rμ, total, seeded, validated}
"""
function cover_rectangle_grid_rect(κ_lo::Float64, κ_hi::Float64,
    μ_lo::Float64, μ_hi::Float64;
    nx::Int, ny::Int,
    bracket_n::Int=400,
    verbose::Bool=true)

    @assert κ_hi > κ_lo "κ_hi must be > κ_lo"
    @assert μ_hi > μ_lo "μ_hi must be > μ_lo"
    @assert nx > 0 && ny > 0 "nx, ny must be positive"

    Lκ = κ_hi - κ_lo
    Lμ = μ_hi - μ_lo
    Δκ = Lκ / nx
    Δμ = Lμ / ny
    rκ = Δκ / 2
    rμ = Δμ / 2

    boxes = Tuple{Bool,Float64,Float64,Float64,Float64}[]
    total = nx * ny
    seeded = 0
    validated = 0

    # iterate cell centers
    for i in 1:nx, j in 1:ny
        κc = κ_lo + (i - 0.5) * Δκ
        μc = μ_lo + (j - 0.5) * Δμ

        # 1) try to bracket a root at this center
        br = brackets_f₁(κc, μc; n=bracket_n)
        if isempty(br)
            continue
        end
        seeded += 1

        # 2) refine the root (Brent on θ)
        root = find_zero(θ -> f₁(θ, κc, μc), br[1])

        # 3) validate this fixed box (no growth)
        ok = validate_f1_root(root, κc, μc, rκ, rμ; verbose=false)
    
        keep_invalid = true  # or make it a keyword
        if ok
            push!(boxes, (true, κc, μc, rκ, rμ))
            validated += 1
        elseif keep_invalid
            push!(boxes, (false, κc, μc, rκ, rμ))
        end

    end

    if verbose
        println("cover_rectangle_grid_rect:")
        println("  Δκ=$(Δκ), Δμ=$(Δμ); rκ=$(rκ), rμ=$(rμ)")
        println("  cells total=$(total), seeded=$(seeded), validated=$(validated)")
    end

    stats = (Δκ=Δκ, Δμ=Δμ, rκ=rκ, rμ=rμ, total=total, seeded=seeded, validated=validated)
    return boxes, stats
end


# Rect(κlo, κhi, μlo, μhi, ok=false, level=0)
function cover_rectangle_grid_rect(R::Rect; nx::Int, ny::Int, bracket_n::Int=400, verbose::Bool=true)
    cover_rectangle_grid_rect(R.κlo, R.κhi, R.μlo, R.μhi;
        nx=nx, ny=ny, bracket_n=bracket_n, verbose=verbose)
end


κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

# e.g., ~0.001×0.001 inside full rectangle → pick nx,ny accordingly
nx, ny = 2464, 89  # or your preferred counts
boxes, stats = cover_rectangle_grid_rect(cells[1]; nx=nx, ny=ny, verbose=true)

println("Validated boxes: ", length(boxes))
# plot_validated_boxes(boxes)
# only invalid boxes
plot_validated_boxes(boxes; mode=:invalid, color_invalid=:red, title="Invalid boxes")


plot_validated_boxes(boxes; mode=:invalid, y_exaggerate=30, plot_centers_invalid=true)

using Printf
# Select invalid boxes from boxes :: Vector{Tuple{Bool,Float64,Float64,Float64,Float64}}
invalid_boxes(boxes) = filter(b -> !b[1], boxes)

# Quick print summary (and peek at the first few)
function show_invalid_summary(boxes; preview=10)
    bad = invalid_boxes(boxes)
    println("Invalid boxes: ", length(bad))
    for (k, b) in enumerate(bad[1:min(preview, end)])
        ok, κ, μ, rκ, rμ = b
        @printf "%3d) κ=[%.6g, %.6g], μ=[%.6g, %.6g]\n" k  (κ-rκ) (κ+rκ) (μ-rμ) (μ+rμ)
    end
    return bad
end


function plot_invalid_centers(boxes; ms=3, color=:red, title="Invalid boxes (centers)")
    bad = invalid_boxes(boxes)
    κs = getindex.(bad, 2)   # centers
    μs = getindex.(bad, 3)
    μT = 0.0886765
    κ_lo, κ_hi = 0.0, (0.5 * π)^2
    μ_lo, μ_hi = 0.0, μT
    Plots.scatter(κs, μs; ms=ms, color=color, xlabel="κ", ylabel="μ",
            title=title, legend=false, xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))
end

bad = show_invalid_summary(boxes)
plot_invalid_centers(boxes)            # simple + visible
