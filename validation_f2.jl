using IntervalArithmetic
using Plots
using Roots
using Printf

include("functions.jl")
include("functions_f2.jl")
# include("validation_functions.jl")
include("plotting_functions.jl")

function brackets_f₂θ(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_2(κ, μ)[1]
    θmax = Intθ_2(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₂θ(θ, κ, μ), a, b; n=n)
end


function classify_f2_rootcount(IθM::Interval, Iκ::Interval, Iμ::Interval;
                               verbose::Bool=false)

    # --- Step 2: θ_M must lie in [√κ q₂, π/2] for all (κ,μ) in the box ---

    Iq2   = q₂(Iκ, Iμ)
    Isqrt = sqrt(Iκ)
    Ileft = Isqrt * Iq2                # all possible √κ q₂

    Iright = interval(pi / 2)          # thin interval {π/2}

    # rigorous allowed range: [sup(Ileft), π/2]
    lower = sup(Ileft)
    upper = sup(Iright)                # = π/2

    # if the lower bound already exceeds π/2, there is no allowed region
    if lower > upper
        verbose && println("No allowed θ-range: [√κ q₂, π/2] is empty.")
        If2_dummy = interval(0)
        return true, :zero_roots, IθM, If2_dummy, interval(lower, upper)
    end

    Iallowed = interval(lower, upper)

    inside = (inf(IθM) ≥ lower) && (sup(IθM) ≤ upper)

    if verbose
        println("θ_M ∈ ", IθM, ", allowed θ ∈ ", Iallowed,
                "  → inside? ", inside)
    end

    if !inside
        # already implies zero roots in this parameter box
        If2_dummy = interval(0)
        return true, :zero_roots, IθM, If2_dummy, Iallowed
    end

    # --- Step 3: sign of f₂ at θ_M ---------------------------------------

    If₂ = f₂(IθM, Iκ, Iμ)
    tg  = sign_tag(If₂)

    if verbose
        println("f₂(θ_M) ∈ ", If₂, "  ⇒ ", tg)
    end

    label =
        tg === :pos ? :two_roots  :
        tg === :neg ? :zero_roots :
                      :unknown

    return true, label, IθM, If₂, Iallowed
end

function validate_f2_root(root, κ, μ, rκ, rμ;
    refine_steps::Int=30,
    verbose::Bool=false)

    # --- 1) parameter intervals ---
    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    # θ-interval [0, θ₁_max] (conservative)
    a = Intθ_2(κ, μ)[1]
    b = Intθ_2(κ, μ)[2]
    Iθ = interval(inf(a), sup(b))

    # guard: root ∈ Iθ
    if !(inf(Iθ) ≤ root ≤ sup(Iθ))
        verbose && println("root=$root not inside θ-interval $Iθ")
        return false, :unknown, Iθ, Interval(0), Interval(0), Interval(0)
    end

    # --- 2) existence by sign change at midpoints ---
    leftθ = (inf(Iθ) + root) / 2
    rightθ = (root + sup(Iθ)) / 2
    fL = f₂θ(interval(leftθ), Iκ, Iμ)
    fR = f₂θ(interval(rightθ), Iκ, Iμ)

    left_neg_right_pos = isneg(fL) && ispos(fR)
    left_pos_right_neg = ispos(fL) && isneg(fR)
    exists = left_neg_right_pos || left_pos_right_neg

    if !exists
        verbose && println("Existence not certified by sign change: fL=$fL, fR=$fR")
        return false, :unknown, Iθ, Interval(0), Interval(0), Interval(0)
    end
    J = interval(leftθ, rightθ)

    # # --- 3) uniqueness: derivative sign + complements ---
    # status, J, dI = certify_deriv_sign_centered(root, κ, μ, rκ, rμ; verbose=false)
    # if status === :unknown
    #     verbose && println("Could not certify strict sign of f₁θ near root. dI = $dI")
    #     return false, :unknown, Iθ, Interval(0), Interval(0), Interval(0)
    # end

    # uniq = validate_complements(Iθ, J, κ, μ, rκ, rμ; verbose=verbose)
    # if !uniq
    #     verbose && println("Uniqueness via complements failed.")
    #     return false, :unknown, J, Interval(0), Interval(0), Interval(0)
    # end

    # verbose && begin
    #     println("θ-interval Iθ     = $Iθ")
    #     println("leftθ=$leftθ  fL  = $fL  (neg? ", isneg(fL), ", pos? ", ispos(fL), ")")
    #     println("rightθ=$rightθ fR = $fR  (neg? ", isneg(fR), ", pos? ", ispos(fR), ")")
    #     println("certified J       = $J (status = $status)")
    #     println("uniqueness via complements? ", uniq)
    # end

    # --- 4) refine the existence interval around root (cheap) ---
    J_ref = refine_root_interval_f2(root, κ, μ, rκ, rμ, J;
        steps=refine_steps, verbose=verbose)

    # --- 5) classify stability on J_ref ---
    okc, label, IθM, If₂, Iallowed = classify_f2_rootcount(J_ref, Iκ, Iμ; verbose=verbose)

    # ok = existence+uniqueness succeeded, classification ran (even if label=:unknown)
    ok_total =  okc

    return ok_total, label, J_ref, IθM, If₂, Iallowed
end

function refine_root_interval_f2(root, κ, μ, rκ, rμ, J::Interval;
    steps::Int=10, verbose::Bool=false)

    Iκ = interval(κ, rκ; format=:midpoint)
    Iμ = interval(μ, rμ; format=:midpoint)

    left = inf(J)
    right = sup(J)

    # ensure root is inside J (it should be, if earlier code is correct)
    if !(left ≤ root ≤ right)
        verbose && println("root=$root not inside certified J=$J; no refinement.")
        return J
    end

    for k in 1:steps
        width = right - left
        width ≤ 0 && break

        # propose a symmetric shrink around root
        δ = width / 4          # shrink factor; can tune
        new_left = root - δ
        new_right = root + δ

        # keep new interval inside original J
        if new_left <= left || new_right >= right
            verbose && println("Step $k: proposed interval not inside J; stop.")
            break
        end

        fL = f₂θ(interval(new_left), Iκ, Iμ)
        fR = f₂θ(interval(new_right), Iκ, Iμ)

        left_neg_right_pos = isneg(fL) && ispos(fR)
        left_pos_right_neg = ispos(fL) && isneg(fR)
        has_sign_change = left_neg_right_pos || left_pos_right_neg

        if has_sign_change
            verbose && println("Step $k: shrink to [$new_left, $new_right], fL=$fL, fR=$fR")
            left, right = new_left, new_right
        else
            verbose && println("Step $k: cannot shrink further (no clean sign change).")
            break
        end
    end

    Jref = interval(left, right)
    verbose && println("Refined J: ", Jref)
    return Jref
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
        br = brackets_f₂θ(κc, μc; n=bracket_n)
        if isempty(br)
            if keep_invalid
                push!(boxes, (false, κc, μc, rκ, rμ, :unknown))
            end
            continue
        end
        seeded += 1

        # 2) refine the root (Brent on θ)
        root = find_zero(θ -> f₂θ(θ, κc, μc), br[1])

        # 3) validate + stability (existence, uniqueness, refined J, classification)
        ok, label, J_ref, IθM, If₂, Iallowed =
            validate_f2_root(root, κc, μc, rκ, rμ;
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

# regions entries: (κ, μ, θM, f2_at_θM, label)

function rootcount_types(regions)
    labels = getindex.(regions, 6)  # grab the 5th element (label)
    return unique(labels)
end

"""
    plot_rootcount_centers(regions;
        μT=0.0886765,
        zoom=nothing,          # (κ_lo, κ_hi, μ_lo, μ_hi) or nothing
        markersize=4,
        title="f₂ root count (centers)",
        colors=Dict(
            :two_roots       => :green,
            :zero_roots      => :red,
            :outside_range   => :gray,
            :no_max          => :black,
            :unknown         => :orange
        ),
        labels=nothing)

Plot the (κ, μ) points, colored by the root-count label from `regions`.

`regions` entries are `(κ::Float64, μ::Float64, θM::Float64, f2val::Float64, label::Symbol)`.

- `labels = nothing`     → show all labels
- `labels = :two_roots`  → only that label
- `labels = [:two_roots, :zero_roots]` → only those labels
"""
function plot_rootcount_centers(regions;
    μT=0.0886765,
    zoom=nothing,
    markersize=4,
    title::AbstractString="f₂ root count (centers)",
    colors=Dict(
        :two_roots       => :green,
        :zero_roots      => :red,
        :outside_range   => :gray,
        :no_max          => :black,
        :unknown         => :orange
    ),
    labels=nothing)

    # normalize label filter to Set{Symbol} or nothing
    labelset = if labels === nothing
        nothing
    elseif labels isa Symbol
        Set{Symbol}((labels,))
    elseif labels isa AbstractVector
        Set(Symbol.(labels))
    else
        error("labels must be nothing, a Symbol, or a vector of Symbols.")
    end

    # full default domain
    κ_lo_full, κ_hi_full = 0.0, (0.5 * π)^2
    μ_lo_full, μ_hi_full = 0.0, μT

    if zoom === nothing
        κ_lo, κ_hi, μ_lo, μ_hi = κ_lo_full, κ_hi_full, μ_lo_full, μ_hi_full
    else
        κ_lo, κ_hi, μ_lo, μ_hi = zoom
    end

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        title=title, legend=true,
        xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    # group by label
    bylabel = Dict{Symbol,Vector{Tuple{Float64,Float64}}}()

    for (κ, μ, θM, f2val, label) in regions
        if labelset !== nothing && !(label in labelset)
            continue
        end
        push!(get!(bylabel, label, Tuple{Float64,Float64}[]), (κ, μ))
    end

    # plot each label with its color
    for (lbl, pts) in sort(collect(bylabel); by = x -> String(x[1]))
        κs = [p[1] for p in pts]
        μs = [p[2] for p in pts]
        col = get(colors, lbl, :black)
        Plots.scatter!(plt, κs, μs;
            color=col, ms=markersize, alpha=0.9,
            label=String(lbl))
    end

    Plots.display(plt)
    return plt
end




nx, ny = 3000, 200
boxes, stats = cover_rectangle_grid_rect(cells[1]; nx=nx, ny=ny, verbose=true)




rootcount_types(boxes)               

# plot_rootcount_centers(boxes)        # full plot
# plot_rootcount_centers(boxes;        # only two-roots region
#     labels=:two_roots,
#     title="f₂: region with two roots")


# plt = plot_invalid_centers(boxes; title ="Invalid boxes(centers) f_2")
# for i in 2:9
#     boxes, stats = cover_rectangle_grid_rect(cells[i]; nx=nx, ny=ny, verbose=true) 
# plot_invalid_centers!(plt, boxes)
# end
# savefig(plt, "testing_grid_f1_invalid_boxes.png")