using IntervalArithmetic
using Plots
using Roots
using Printf

include("functions.jl")
include("validation_functions.jl")
include("plotting_functions.jl")

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

# Domain defaults: κ ∈ [0, (π/2)^2], μ ∈ [0, μT]
κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

inner, cells = inner_3x3(κ_lo, κ_hi, μ_lo, μ_hi; margin=0.002)
println("Inner rect: ", inner)
println("3×3 cells = ", length(cells))
println("First cell: ", cells[1])
println("Last  cell: ", cells[end])



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

function cover_rectangle_grid_rect(R::Rect; nx::Int, ny::Int, bracket_n::Int=400, verbose::Bool=true)
    cover_rectangle_grid_rect(R.κlo, R.κhi, R.μlo, R.μhi;
        nx=nx, ny=ny, bracket_n=bracket_n, verbose=verbose)
end

κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

# nx, ny = 2464, 89
# boxes, stats = cover_rectangle_grid_rect(cells[1]; nx=nx, ny=ny, verbose=true)

# plt = plot_invalid_centers(boxes; title ="Invalid boxes(centers) f_1")
# for i in 2:9
#     boxes, stats = cover_rectangle_grid_rect(cells[i]; nx=nx, ny=ny, verbose=true) 
# plot_invalid_centers!(plt, boxes)
# end
# savefig(plt, "testing_grid_f1_invalid_boxes.png")

### Stability functions
function cover_rectangle_grid_rect(κ_lo::Float64, κ_hi::Float64,
    μ_lo::Float64, μ_hi::Float64;
    nx::Int, ny::Int,
    bracket_n::Int=400,
    refine_steps::Int=8,
    keep_invalid::Bool=true,
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

    # now includes stability label :: Symbol
    boxes = Tuple{Bool,Float64,Float64,Float64,Float64,Symbol}[]
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
            if keep_invalid
                push!(boxes, (false, κc, μc, rκ, rμ, :unknown))
            end
            continue
        end
        seeded += 1

        # 2) refine the root (Brent on θ)
        root = find_zero(θ -> f₁(θ, κc, μc), br[1];atol=1e-12)

        # 3) validate + stability (existence, uniqueness, refined J, classification)
        ok, label, J_ref, Idet, Ilam2, IΔ =
            validate_f1_root_with_stability(root, κc, μc, rκ, rμ;
                refine_steps=refine_steps,
                verbose=false)

        if ok
            push!(boxes, (true, κc, μc, rκ, rμ, label))
            validated += 1
        elseif keep_invalid
            push!(boxes, (false, κc, μc, rκ, rμ, label))
        end
    end

    if verbose
        println("cover_rectangle_grid_rect (with stability):")
        println("  Δκ=$(Δκ), Δμ=$(Δμ); rκ=$(rκ), rμ=$(rμ)")
        println("  cells total=$(total), seeded=$(seeded), validated=$(validated)")
    end

    stats = (Δκ=Δκ, Δμ=Δμ, rκ=rκ, rμ=rμ, total=total, seeded=seeded, validated=validated)
    return boxes, stats
end

function cover_rectangle_grid_rect(R::Rect;
    nx::Int, ny::Int,
    bracket_n::Int=400,
    refine_steps::Int=8,
    keep_invalid::Bool=true,
    verbose::Bool=true)

    cover_rectangle_grid_rect(R.κlo, R.κhi, R.μlo, R.μhi;
        nx=nx, ny=ny,
        bracket_n=bracket_n,
        refine_steps=refine_steps,
        keep_invalid=keep_invalid,
        verbose=verbose)
end

nx, ny = 2464, 89
boxes, stats = cover_rectangle_grid_rect(cells[9]; nx=nx, ny=ny, verbose=true)
plt = plot_stability_centers(boxes; labels=:elliptic)

# for i in 2:9
#     boxes, stats = cover_rectangle_grid_rect(cells[i]; nx=nx, ny=ny, verbose=true) 
# plot_stability_centers!(plt, boxes; labels=:elliptic)
# end

# savefig(plt, "testing_grid_f1_stability_centers.png")


stability_types(boxes)