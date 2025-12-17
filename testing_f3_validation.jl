using IntervalArithmetic
using Plots
# using PlotlyJS
using Roots
using Printf

include("functions.jl")
include("stability_functions.jl")
include("functions_f3.jl")
include("validation_functions2.jl")
include("plotting_functions.jl")


function classify_f3_rootcount(Iθm::Interval, Iκ::Interval, Iμ::Interval; verbose::Bool=false)
    # --- Step 2: θ_m must be inside [ -π/2, -√κ q₁ ] ----------------------

    Ileft = interval(-pi / 2)                 # thin interval lower bound
    Iright = -(sqrt(Iκ) * q₁(Iκ, Iμ))         # interval upper bound
    Iallowed = interval(inf(Ileft), sup(Iright))

    inside = (inf(Iθm) ≥ inf(Iallowed)) && (sup(Iθm) ≤ sup(Iallowed))

    if verbose
        println("θ_m ∈ ", Iθm, " ; allowed ∈ ", Iallowed, " ; inside? ", inside)
    end

    if !inside
        # certified outside ⇒ zero roots
        return true, :zero_roots, Iθm, interval(0), Iallowed
    end

    # --- Step 3: sign of f₃ at θ_m ---------------------------------------

    If3 = f₃(Iθm, Iκ, Iμ)
    tg = sign_tag(If3)

    if verbose
        println("f₃(θ_m) ∈ ", If3, " ⇒ ", tg)
    end

    if tg === :neg
        return true, :two_roots, Iθm, If3, Iallowed
    elseif tg === :pos
        return true, :zero_roots, Iθm, If3, Iallowed
    else
        return false, :unknown, Iθm, If3, Iallowed
    end
end

function θm_from_br(br, κ, μ; atol=1e-12, verbose=false)
    bestθ = NaN
    bestval = Inf
    found = false

    for x0 in br
        θ = find_zero(θ -> f₃θ(θ, κ, μ), x0; atol=atol)

        s2 = f₃θθ(θ, κ, μ)      # second derivative test
        if s2 > 0               # local minimum
            val = f₃(θ, κ, μ)   # height of the minimum
            if !found || val < bestval
                bestval = val
                bestθ = θ
                found = true
            end
        end

        verbose && println("x0=$x0  θ*=$θ  f₃θθ=$s2  f₃=$(f₃(θ, κ, μ))")
    end

    return found ? bestθ : NaN
end


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

    boxes = Tuple{Bool,Float64,Float64,Float64,Float64,Symbol}[]
    total = nx * ny
    seeded = 0
    validated = 0

    # iterate cell centers
    for i in 1:nx, j in 1:ny
        κc = κ_lo + (i - 0.5) * Δκ
        μc = μ_lo + (j - 0.5) * Δμ

        (θmin, θmax) = Intθ_3(κc, μc)

        br = brackets_signchanges(θ -> f₃θ(θ, κc, μc), θmin, θmax)

        if isempty(br)
            continue
        end
        seeded += 1

        # 2) refine the root (Brent on θ)
        root = θm_from_br(br, κc, μc; atol=1e-12)


        Iκ = interval(κc, rκ; format=:midpoint)
        Iμ = interval(μc, rμ; format=:midpoint)

        ok, data = validate_root_sign_mono(
            f₃θ,
            f₃θθ,
            root,
            interval(θmin, θmax),
            Iκ,
            Iμ;
            ntiles=16,
            min_width=1e-12,
            min_radius=1e-12,
            verbose=false,
        )

        keep_invalid = true  # or make it a keyword
        if ok
            clasif = classify_f3_rootcount(data.Iexist, Iκ, Iμ; verbose=false)
            push!(boxes, (true, κc, μc, rκ, rμ, clasif[2]))
            validated += 1
        elseif keep_invalid
            push!(boxes, (false, κc, μc, rκ, rμ, :unknown))
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
nx, ny = 3000, 100
boxes, stats = cover_rectangle_grid_rect(inner; nx=nx, ny=ny, verbose=true)

using Plots

# --- modern palette (two colors + gray) ---
@inline function rootcount_color(lbl::Symbol)
    lbl === :two_roots && return "#3B6FB6"  # deep blue
    lbl === :zero_roots && return "#E69F00"  # warm amber
    lbl === :unknown && return "#ED1B2F"
    return "#000000"
end

"""
Allocating: new plot of centers colored by label.
boxes entries: (ok, κc, μc, rκ, rμ, label::Symbol)
"""
function plot_rootcount_centers_alloc(boxes;
    only_ok::Bool=true,
    labels=nothing,
    ms::Real=1.0,
    alpha::Real=0.15,
    title::AbstractString="Root-count labels",
    xlabel::AbstractString="κ",
    ylabel::AbstractString="μ",
    legend::Bool=false,
    framestyle::Symbol=:box,
    spinewidth::Real=2.5,
    tickwidth::Real=2.0,
)

    labelset = labels === nothing ? nothing :
               labels isa Symbol ? Set{Symbol}((labels,)) :
               Set(Symbol.(labels))

    κvals = Float64[]
    μvals = Float64[]
    cols = String[]   # store hex strings

    for (ok, κc, μc, rκ, rμ, lbl) in boxes
        if only_ok && !ok
            continue
        end
        if labelset !== nothing && !(lbl in labelset)
            continue
        end
        push!(κvals, κc)
        push!(μvals, μc)
        push!(cols, rootcount_color(lbl))
    end

    plt = Plots.scatter(κvals, μvals;
        ms=ms, color=cols,
        alpha=alpha,
        xlabel=xlabel, ylabel=ylabel,
        title=title, legend=legend,
        framestyle=framestyle,
        foreground_color_border=:black,
        foreground_color_axis=:black,
        foreground_color_text=:black,
        linewidth=spinewidth,
        ticklinewidth=tickwidth,
        grid=false,
    )

    return plt
end

# 
plt = plot_rootcount_centers_alloc(boxes; title="Root-count labels (centers) f_3")
# for i in 1:9
#     boxes, stats = cover_rectangle_grid_rect(cells[i]; nx=nx, ny=ny, verbose=true)
#     plot_rootcount_centers!(plt, boxes)
# end
savefig(plt, "testing_grid_f3_invalid_boxes.png")


