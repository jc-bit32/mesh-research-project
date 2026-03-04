using Comodo

R = 1.0
subdiv = 0  # very coarse

F, V = Comodo.geosphere(subdiv, R)

nverts = length(V)      # V is a Vector of 3D points
nfaces = length(F)

println("Number of vertices = ", nverts)
println("Number of faces    = ", nfaces)

# helper: convert GeometryBasics.Point{3} -> NTuple{3,Float64}
to3(p) = (Float64(p[1]), Float64(p[2]), Float64(p[3]))

# Compute per-face triangle areas
face_area = zeros(nfaces)

for i in 1:nfaces
    i1, i2, i3 = Tuple(F[i])

    x1 = to3(V[i1])
    x2 = to3(V[i2])
    x3 = to3(V[i3])

    # edge vectors
    a = (x2[1]-x1[1], x2[2]-x1[2], x2[3]-x1[3])
    b = (x3[1]-x1[1], x3[2]-x1[2], x3[3]-x1[3])

    # cross(a,b) by hand (so we don't need LinearAlgebra)
    cx = a[2]*b[3] - a[3]*b[2]
    cy = a[3]*b[1] - a[1]*b[3]
    cz = a[1]*b[2] - a[2]*b[1]

    face_area[i] = 0.5 * sqrt(cx*cx + cy*cy + cz*cz)
end

# simple stats
minA = minimum(face_area)
maxA = maximum(face_area)
meanA = sum(face_area) / length(face_area)
totalA = sum(face_area)

println("Face area min  = ", minA)
println("Face area max  = ", maxA)
println("Face area mean = ", meanA)
println("Approx total area = ", totalA)
println("4π = ", 4π)

println("\nFirst few vertices (as x,y,z):")
for i in 1:min(nverts, 20)
    x = to3(V[i])
    println(i, ": ", x)
end