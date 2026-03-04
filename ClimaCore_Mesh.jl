using ClimaCore
using ClimaComms
using Printf


R  = 1.0
ne = 1          # 1 element per face
Nq = 2          # GLL(2) corners



domain = ClimaCore.Domains.SphereDomain(R)
mesh   = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne)

context = ClimaComms.SingletonCommsContext()
topo = ClimaCore.Topologies.Topology2D(context, mesh)

quad  = ClimaCore.Quadratures.GLL{Nq}()
space = ClimaCore.Spaces.SpectralElementSpace2D(topo, quad)

geom = ClimaCore.Fields.local_geometry_field(space)
coords_field = geom.coordinates
J_field      = geom.J


A = Array(parent(coords_field))
@printf("coords parent eltype = %s\n", string(eltype(A)))
@printf("coords parent size   = %s\n", string(size(A)))


lat = vec(@view A[:, :, 1, :])   
lon = vec(@view A[:, :, 2, :])

n = length(lat)


x = R .* cos.(lat) .* cos.(lon)
y = R .* cos.(lat) .* sin.(lon)
z = R .* sin.(lat)

coords = hcat(x, y, z)


J_vec = vec(Array(parent(J_field)))

@printf("npoints (coords) = %d\n", n)
@printf("npoints (J)      = %d\n", length(J_vec))

nuse = min(n, length(J_vec))


println("index, x, y, z, lat(rad), lon(rad), J")
for i in 1:min(nuse, 30)
    @printf("%3d, % .6e, % .6e, % .6e, % .6e, % .6e, % .6e\n",
        i, coords[i,1], coords[i,2], coords[i,3], lat[i], lon[i], J_vec[i])
end

