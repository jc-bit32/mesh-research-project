using ClimaCore
using ClimaComms
using Comodo
using GLMakie
using Rotations
using LinearAlgebra

# idealized 

function equidistant_cubed_sphere_mapping(x, y)
    zx = tan(π/4 * x)
    zy = tan(π/4 * y)
    z0 = 1 / sqrt(1 + zx^2 + zy^2)
    X = zx * z0
    Y = zy * z0
    Z = z0
    return X, Y, Z
end

N = 17
xs = range(-1, 1, length=N)
ys = range(-1, 1, length=N)

X = zeros(Float64, N, N)
Y = zeros(Float64, N, N)
Z = zeros(Float64, N, N)

for j in 1:N, i in 1:N
    X[i, j], Y[i, j], Z[i, j] = equidistant_cubed_sphere_mapping(xs[i], ys[j])
end

rotations = (
    I,
    RotY(-π/2),
    RotY(π/2),
    RotX(-π/2),
    RotX(π/2),
    RotX(π)
)

# ClimaCore

R  = 1.0
ne = 2
Nq = 2

domain = ClimaCore.Domains.SphereDomain(R)
mesh   = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne)

context = ClimaComms.SingletonCommsContext()
topo    = ClimaCore.Topologies.Topology2D(context, mesh)

quad  = ClimaCore.Quadratures.GLL{Nq}()
space = ClimaCore.Spaces.SpectralElementSpace2D(topo, quad)

geom = ClimaCore.Fields.local_geometry_field(space)
A = Array(parent(geom.coordinates))


#Comodo

r = 1.0
nsub = 1

F, V = Comodo.platonicsolid(2, r)
Fn, Vn = Comodo.subquad(F, V, nsub; method = :Catmull_Clark)

for i in 1:length(Vn)
    p = Vn[i]
    s = sqrt(p[1]^2 + p[2]^2 + p[3]^2)
    Vn[i] = p / s
end

vx = [p[1] for p in Vn]
vy = [p[2] for p in Vn]
vz = [p[3] for p in Vn]

# Plotting

fig = Figure(size = (1500, 500))

ax1 = Axis3(fig[1, 1], aspect = :data, title = "Idealized Mesh")
ax2 = Axis3(fig[1, 2], aspect = :data, title = "ClimaCore Coarse Mesh")
ax3 = Axis3(fig[1, 3], aspect = :data, title = "Comodo Coarse Mesh")

# --- idealized wireframe ---
for rot in rotations
    Xp = similar(X)
    Yp = similar(Y)
    Zp = similar(Z)

    for idx in CartesianIndices(X)
    Xp[idx], Yp[idx], Zp[idx] = rot * [X[idx], Y[idx], Z[idx]]
    end

    wireframe!(ax1, Xp, Yp, Zp)
end

# --- ClimaCore full mesh background ---
for e in axes(A, 4)
    lat_deg = @view A[:, :, 1, e]
    lon_deg = @view A[:, :, 2, e]

    lat = deg2rad.(lat_deg)
    lon = deg2rad.(lon_deg)

    Xc = cos.(lat) .* cos.(lon)
    Yc = cos.(lat) .* sin.(lon)
    Zc = sin.(lat)

    ni, nj = size(Xc)

    lines!(ax2, Xc[:, 1],  Yc[:, 1],  Zc[:, 1], color = (:gray, 0.35), linewidth = 1)
    lines!(ax2, Xc[:, nj], Yc[:, nj], Zc[:, nj], color = (:gray, 0.35), linewidth = 1)
    lines!(ax2, Xc[1, :],  Yc[1, :],  Zc[1, :], color = (:gray, 0.35), linewidth = 1)
    lines!(ax2, Xc[ni, :], Yc[ni, :], Zc[ni, :], color = (:gray, 0.35), linewidth = 1)
end

# --- raw ClimaCore storage order, first occurrence only ---
lat_all_deg = vec(@view A[:, :, 1, :])
lon_all_deg = vec(@view A[:, :, 2, :])

lat_all = deg2rad.(lat_all_deg)
lon_all = deg2rad.(lon_all_deg)

xc_all = cos.(lat_all) .* cos.(lon_all)
yc_all = cos.(lat_all) .* sin.(lon_all)
zc_all = sin.(lat_all)

seen = Set{NTuple{3,Int}}()
keep_idx = Int[]

for k in eachindex(xc_all)
    key = (
        round(Int, xc_all[k] * 10^8),
        round(Int, yc_all[k] * 10^8),
        round(Int, zc_all[k] * 10^8),
    )
    if !(key in seen)
        push!(seen, key)
        push!(keep_idx, k)
    end
end

xk = xc_all[keep_idx]
yk = yc_all[keep_idx]
zk = zc_all[keep_idx]

# lighter path + readable labels
scatter!(ax2, xk, yk, zk, markersize = 10, color = :dodgerblue)
lines!(ax2, xk, yk, zk, color = :darkorange, linewidth = 2)

for n in eachindex(keep_idx)
    text!(
        ax2,
        xk[n] * 1.03, yk[n] * 1.03, zk[n] * 1.03,
        text = string(n),
        fontsize = 10
    )
end

# --- comodo points + labels ---
scatter!(ax3, vx, vy, vz, markersize = 16)

for i in eachindex(vx)
    text!(ax3, vx[i], vy[i], vz[i], text = string(i), fontsize = 14)
end

# --- comodo quad edges ---
for face in Fn
    i1, i2, i3, i4 = Tuple(face)

    p1 = Vn[i1]
    p2 = Vn[i2]
    p3 = Vn[i3]
    p4 = Vn[i4]

    lines!(ax3, [p1[1], p2[1]], [p1[2], p2[2]], [p1[3], p2[3]])
    lines!(ax3, [p2[1], p3[1]], [p2[2], p3[2]], [p2[3], p3[3]])
    lines!(ax3, [p3[1], p4[1]], [p3[2], p4[2]], [p3[3], p4[3]])
    lines!(ax3, [p4[1], p1[1]], [p4[2], p1[2]], [p4[3], p1[3]])
end

display(fig)
save("combined_mesh_debug.png", fig)