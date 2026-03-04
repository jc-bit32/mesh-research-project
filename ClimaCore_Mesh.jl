using ClimaCore
using ClimaComms

R  = 1.0
ne = 1
Nq = 2

domain = ClimaCore.Domains.SphereDomain(R)
mesh   = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne)

context = ClimaComms.SingletonCommsContext()
topo    = ClimaCore.Topologies.Topology2D(context, mesh)

quad  = ClimaCore.Quadratures.GLL{Nq}()
space = ClimaCore.Spaces.SpectralElementSpace2D(topo, quad)

geom = ClimaCore.Fields.local_geometry_field(space)

A = Array(parent(geom.coordinates))   # (Nq, Nq, 2, 6) in your build
lat = vec(@view A[:, :, 1, :])
lon = vec(@view A[:, :, 2, :])

x = R .* cos.(lat) .* cos.(lon)
y = R .* cos.(lat) .* sin.(lon)
z = R .* sin.(lat)

J = vec(Array(parent(geom.J)))

n = min(length(J), length(lat))

println("index,x,y,z,lat,lon,J")
for i in 1:n
    println(i, ",", x[i], ",", y[i], ",", z[i], ",", lat[i], ",", lon[i], ",", J[i])
end