using Comodo

r = 1.0
nsub = 1

# cube -> Catmull-Clark subdivision
F, V = Comodo.platonicsolid(2, r)
Fn, Vn = Comodo.subquad(F, V, nsub; method = :Catmull_Clark)

# normalize to unit sphere
for i in 1:length(Vn)
    p = Vn[i]
    s = sqrt(p[1]^2 + p[2]^2 + p[3]^2)
    Vn[i] = p / s
end

println("verts = ", length(Vn))
println("faces = ", length(Fn))

# radius sanity check
radii = [sqrt(p[1]^2 + p[2]^2 + p[3]^2) for p in Vn]
println("radius min/max = ", minimum(radii), " / ", maximum(radii))

println("\nfirst few vertices:")
for i in 1:min(12, length(Vn))
    p = Vn[i]
    println(i, ": (", p[1], ", ", p[2], ", ", p[3], ")")
end