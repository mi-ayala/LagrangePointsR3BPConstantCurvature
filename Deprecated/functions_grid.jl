# make an interior grid (with margin)
interior_linspace(lo, hi, n; margin=0.02) = begin
    a = lo + margin * (hi - lo)
    b = hi - margin * (hi - lo)
    range(a, b; length=n)
end