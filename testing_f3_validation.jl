using IntervalArithmetic
using Plots
# using PlotlyJS
using Roots
using Printf

include("functions_helpers.jl")
include("functions_stability.jl")
include("functions_f3.jl")
include("functions_validation.jl")
include("functions_plots.jl")


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
# plt = plot_rootcount_centers_alloc(boxes; title="Root-count labels (centers) f_3")
# for i in 1:9
#     boxes, stats = cover_rectangle_grid_rect(cells[i]; nx=nx, ny=ny, verbose=true)
#     plot_rootcount_centers!(plt, boxes)
# end
# savefig(plt, "testing_grid_f3_invalid_boxes.png")

function rasterize_boxes_by_label(boxes;
    want::Symbol,
    κ_lo=0.0, κ_hi=(0.5*pi)^2,
    μ_lo=0.0, μ_hi=0.0886765,
    nx::Int=800, ny::Int=800,
    only_ok::Bool=true)

    κgrid = collect(range(κ_lo, κ_hi; length=nx))
    μgrid = collect(range(μ_lo, μ_hi; length=ny))
    occ = falses(nx, ny)

    # quick helpers: map coordinate → index
    @inline function ix(x)
        t = (x - κ_lo) / (κ_hi - κ_lo)
        return clamp(Int(floor(t*(nx-1))) + 1, 1, nx)
    end
    @inline function iy(y)
        t = (y - μ_lo) / (μ_hi - μ_lo)
        return clamp(Int(floor(t*(ny-1))) + 1, 1, ny)
    end

    for (ok, κc, μc, rκ, rμ, lbl) in boxes
        only_ok && !ok && continue
        lbl === want || continue

        κ1 = max(κ_lo, κc - rκ); κ2 = min(κ_hi, κc + rκ)
        μ1 = max(μ_lo, μc - rμ); μ2 = min(μ_hi, μc + rμ)

        i1, i2 = ix(κ1), ix(κ2)
        j1, j2 = iy(μ1), iy(μ2)

        @inbounds occ[i1:i2, j1:j2] .= true
    end

    return κgrid, μgrid, occ
end

function boundary_mask(occ::AbstractMatrix{Bool})
    nx, ny = size(occ)
    b = falses(nx, ny)

    @inbounds for i in 2:nx-1, j in 2:ny-1
        if occ[i,j] && (!occ[i-1,j] || !occ[i+1,j] || !occ[i,j-1] || !occ[i,j+1])
            b[i,j] = true
        end
    end

    # edges: treat any true edge cell as boundary
    @inbounds for i in 1:nx
        b[i,1]  |= occ[i,1]
        b[i,ny] |= occ[i,ny]
    end
    @inbounds for j in 1:ny
        b[1,j]  |= occ[1,j]
        b[nx,j] |= occ[nx,j]
    end

    return b
end

function plot_rootcount_boundaries(boxes;
    κ_lo=0.0, κ_hi=(0.5*pi)^2,
    μ_lo=0.0, μ_hi=0.0886765,
    nx::Int=800, ny::Int=800,
    only_ok::Bool=true,
    ms::Real=0.6,
    alpha::Real=0.25)

    plt = Plots.plot(; xlabel="κ", ylabel="μ", legend=true,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi),
        framestyle=:box)

    # change these labels to match YOUR symbols
    for (lbl, col) in ((:one_root, "#2D6A4F"), (:two_roots, "#3A86FF"))
        κgrid, μgrid, occ = rasterize_boxes_by_label(
            boxes; want=lbl,
            κ_lo=κ_lo, κ_hi=κ_hi, μ_lo=μ_lo, μ_hi=μ_hi,
            nx=nx, ny=ny, only_ok=only_ok
        )
        b = boundary_mask(occ)
        I = findall(b)

        κb = [κgrid[t[1]] for t in I]
        μb = [μgrid[t[2]] for t in I]

        Plots.scatter!(plt, κb, μb; ms=ms, alpha=alpha, color=col, label=String(lbl))
    end

    return plt
end

plt = plot_rootcount_boundaries(boxes; nx=1000, ny=1000, alpha=0.25, ms=0.5)
Plots.display(plt)

function plot_rootcount_fill_two(boxes;
    κ_lo=0.0, κ_hi=(0.5*pi)^2,
    μ_lo=0.0, μ_hi=0.0886765,
    nx::Int=600, ny::Int=600,
    only_ok::Bool=true,
    lbl1::Symbol=:two_roots, col1="#3A86FF",
    lbl2::Symbol=:zero_roots, col2="#2D6A4F",
    alpha1::Real=0.25, alpha2::Real=0.20)

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi),
        framestyle=:box, legend=false)

    # first fill
    κg, μg, occ = rasterize_boxes_by_label(boxes; want=lbl1,
        κ_lo=κ_lo, κ_hi=κ_hi, μ_lo=μ_lo, μ_hi=μ_hi,
        nx=nx, ny=ny, only_ok=only_ok)
    Z = Float64.(occ)
    Plots.heatmap!(plt, κg, μg, Z'; color=[RGBA(0,0,0,0), Plots.parse(Colorant, col1)],
        clims=(0.0, 1.0), colorbar=false, alpha=alpha1)

    # second fill
    κg, μg, occ = rasterize_boxes_by_label(boxes; want=lbl2,
        κ_lo=κ_lo, κ_hi=κ_hi, μ_lo=μ_lo, μ_hi=μ_hi,
        nx=nx, ny=ny, only_ok=only_ok)
    Z = Float64.(occ)
    Plots.heatmap!(plt, κg, μg, Z'; color=[RGBA(0,0,0,0), Plots.parse(Colorant, col2)],
        clims=(0.0, 1.0), colorbar=false, alpha=alpha2)

    return plt
end

plt = plot_rootcount_fill_two(boxes;
    lbl1=:two_roots, col1="#3A86FF",
    lbl2=:zero_roots, col2="#2D6A4F",
    alpha1=0.7, alpha2=0.7,
    nx=900, ny=900)
Plots.display(plt)

function plot_rootcount_fill_three(boxes;
    κ_lo=0.0, κ_hi=(0.5*pi)^2,
    μ_lo=0.0, μ_hi=0.0886765,
    nx::Int=600, ny::Int=600,
    only_ok::Bool=true,
    lbl1::Symbol=:two_roots,   col1="#3A86FF", alpha1::Real=0.4,
    lbl2::Symbol=:zero_roots,  col2="#2D6A4F", alpha2::Real=0.4,
    lbl3::Symbol=:unknown,     col3="#ED1B2F", alpha3::Real=0.4,
    title::AbstractString="Root-count fill")

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi),
        framestyle=:box, legend=false, title=title)

    # helper to overlay one label
    function add_layer!(lbl, col, α)
        κg, μg, occ = rasterize_boxes_by_label(
            boxes; want=lbl,
            κ_lo=κ_lo, κ_hi=κ_hi, μ_lo=μ_lo, μ_hi=μ_hi,
            nx=nx, ny=ny, only_ok=only_ok
        )
        Z = Float64.(occ)
        Plots.heatmap!(plt, κg, μg, Z';
            color=[RGBA(0,0,0,0), Plots.parse(Colorant, col)],
            clims=(0.0, 1.0),
            colorbar=false,
            alpha=α)
        return nothing
    end

    add_layer!(lbl1, col1, alpha1)
    add_layer!(lbl2, col2, alpha2)
    add_layer!(lbl3, col3, alpha3)

    return plt
end

plt = plot_rootcount_fill_three(boxes; nx=600, ny=600, col3="#ED1B2F", lbl3=:unknown)
Plots.display(plt)