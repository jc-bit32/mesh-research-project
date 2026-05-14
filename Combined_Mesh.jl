using ClimaCore
using ClimaComms
using Comodo
using Statistics
using CairoMakie
using LinearAlgebra
using Revise

angle_diff_deg(a, b) = abs(mod(a - b + 180, 360) - 180)

function compare_latlon_fields(native_field, test_field)
    A1 = Array(parent(native_field))
    A2 = Array(parent(test_field))

    max_err = 0.0
    worst = (0, 0, 0, "none")

    nelems = size(A1, 4)

    for e in 1:nelems
        for j in axes(A1, 2), i in axes(A1, 1)
            lat_err = abs(A1[i, j, 1, e] - A2[i, j, 1, e])
            lon_err = angle_diff_deg(A1[i, j, 2, e], A2[i, j, 2, e])

            if lat_err > max_err
                max_err = lat_err
                worst = (i, j, e, "lat")
            end
            if lon_err > max_err
                max_err = lon_err
                worst = (i, j, e, "lon")
            end
        end
    end

    return max_err, worst
end

function compare_numeric_fields(native_field, test_field)
    A1 = Array(parent(native_field))
    A2 = Array(parent(test_field))
    diff = abs.(A1 .- A2)
    return maximum(diff), size(A1)
end

function stats_summary(v)
    s = sort(collect(v))
    return (
        min    = minimum(s),
        q1     = quantile(s, 0.25),
        median = median(s),
        mean   = mean(s),
        q3     = quantile(s, 0.75),
        max    = maximum(s),
        std    = std(s),
    )
end

function print_stats(name, s)
    println(name)
    println("  min    = ", s.min)
    println("  q1     = ", s.q1)
    println("  median = ", s.median)
    println("  mean   = ", s.mean)
    println("  q3     = ", s.q3)
    println("  max    = ", s.max)
    println("  std    = ", s.std)
end

function clima_element_corners(clima_xyz, e, Nq)
    return [
        vec(clima_xyz[1,  1,  :, e]),
        vec(clima_xyz[Nq, 1,  :, e]),
        vec(clima_xyz[Nq, Nq, :, e]),
        vec(clima_xyz[1,  Nq, :, e]),
    ]
end

function comodo_face_corners(Vn, Fn, f)
    i1, i2, i3, i4 = Tuple(Fn[f])
    return [Vn[i1], Vn[i2], Vn[i3], Vn[i4]]
end

function comodo_xyz_at_ξη_from_face(Vn, Fn, f, ξ1, ξ2)
    i1, i2, i3, i4 = Tuple(Fn[f])

    r11 = Vn[i1]
    r21 = Vn[i2]
    r22 = Vn[i3]
    r12 = Vn[i4]

    w11 = 0.25 * (1 - ξ1) * (1 - ξ2)
    w21 = 0.25 * (1 + ξ1) * (1 - ξ2)
    w22 = 0.25 * (1 + ξ1) * (1 + ξ2)
    w12 = 0.25 * (1 - ξ1) * (1 + ξ2)

    return w11 .* r11 .+ w21 .* r21 .+ w22 .* r22 .+ w12 .* r12
end

function latlon_block_to_xyz(A, Nq)
    nelems = size(A, 4)
    clima_xyz = zeros(Float64, Nq, Nq, 3, nelems)

    for e in 1:nelems
        lat_deg = @view A[:, :, 1, e]
        lon_deg = @view A[:, :, 2, e]

        lat = deg2rad.(lat_deg)
        lon = deg2rad.(lon_deg)

        clima_xyz[:, :, 1, e] = cos.(lat) .* cos.(lon)
        clima_xyz[:, :, 2, e] = cos.(lat) .* sin.(lon)
        clima_xyz[:, :, 3, e] = sin.(lat)
    end

    return clima_xyz
end

nsub = 2
R = 1.0
ne = 2^nsub
Nq = 2
mapping_types = [:equiangular, :equidistant, :conformal]


#Comodo

r = 1.0

use_exact_match = (Nq == 2 && nsub == 1)

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

# ============================================================
# Build Comodo mesh
# ============================================================

F, V = Comodo.platonicsolid(2, r)
Fn, Vn = Comodo.subquad(F, V, nsub; method = :Catmull_Clark)

for i in 1:length(Vn)
    p = Vn[i]
    s = sqrt(p[1]^2 + p[2]^2 + p[3]^2)
    Vn[i] = p / s
end


for mapping_type in mapping_types

    # ============================================================
    # ClimaCore setup
    # ============================================================


    domain = ClimaCore.Domains.SphereDomain(R)

    if mapping_type == :equiangular
        mesh = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne)
    elseif mapping_type == :equidistant
        mesh = ClimaCore.Meshes.EquidistantCubedSphere(domain, ne)
    elseif mapping_type == :conformal
        mesh = ClimaCore.Meshes.ConformalCubedSphere(domain, ne)
    else
        error("Unknown mapping_type = $mapping_type")
    end

    context = ClimaComms.SingletonCommsContext()
    topo = ClimaCore.Topologies.Topology2D(context, mesh)
    quad = ClimaCore.Quadratures.GLL{Nq}()
    space = ClimaCore.Spaces.SpectralElementSpace2D(topo, quad)
    geom = ClimaCore.Fields.local_geometry_field(space)
    A = Array(parent(geom.coordinates))




    
    # ============================================================
    # Build ClimaCore xyz from native geometry coordinates
    # ============================================================

    A = Array(parent(geom.coordinates))
    nelems = size(A, 4)
    clima_xyz = latlon_block_to_xyz(A, Nq)

    quad_points, quad_weights = ClimaCore.Quadratures.quadrature_points(Float64, quad)
    Dmat = ClimaCore.Quadratures.differentiation_matrix(Float64, quad)

    

    # ============================================================
    # Match each ClimaCore element to a Comodo face using corners
    # ============================================================

    face_match = zeros(Int, nelems)
    used_faces = falses(length(Fn))

    for e in 1:nelems
        cpts = clima_element_corners(clima_xyz, e, Nq)

        best_f = 0
        best_score = Inf

        for f in eachindex(Fn)
            used_faces[f] && continue

            fpts = comodo_face_corners(Vn, Fn, f)

            score = 0.0
            for cp in cpts
                best_d = Inf
                for fp in fpts
                    dx = cp[1] - fp[1]
                    dy = cp[2] - fp[2]
                    dz = cp[3] - fp[3]
                    d2 = dx*dx + dy*dy + dz*dz
                    best_d = min(best_d, d2)
                end
                score += best_d
            end

            if score < best_score
                best_score = score
                best_f = f
            end
        end

        if best_f == 0
            error("No Comodo face match found for ClimaCore element $e")
        end

        face_match[e] = best_f
        used_faces[best_f] = true
    end

    # ============================================================
    # Reconstruct Comodo geometry on ClimaCore nodal locations
    # ============================================================

    comodo_xyz = zeros(Float64, Nq, Nq, 3, nelems)

    if use_exact_match
        

        for e in 1:nelems
            f = face_match[e]
            i1, i2, i3, i4 = Tuple(Fn[f])

            pts = [
                Vn[i1],
                Vn[i2],
                Vn[i3],
                Vn[i4],
            ]

            for j in 1:Nq, i in 1:Nq
                xc = clima_xyz[i, j, 1, e]
                yc = clima_xyz[i, j, 2, e]
                zc = clima_xyz[i, j, 3, e]

                best = 1
                best_d = Inf

                for k in 1:4
                    dx = xc - pts[k][1]
                    dy = yc - pts[k][2]
                    dz = zc - pts[k][3]
                    d2 = dx*dx + dy*dy + dz*dz

                    if d2 < best_d
                        best_d = d2
                        best = k
                    end
                end

                p = pts[best]
                comodo_xyz[i, j, 1, e] = p[1]
                comodo_xyz[i, j, 2, e] = p[2]
                comodo_xyz[i, j, 3, e] = p[3]
            end
        end
    else
        

        for e in 1:nelems
            f = face_match[e]

            for j in 1:Nq, i in 1:Nq
                x = comodo_xyz_at_ξη_from_face(Vn, Fn, f, quad_points[i], quad_points[j])

                comodo_xyz[i, j, 1, e] = x[1]
                comodo_xyz[i, j, 2, e] = x[2]
                comodo_xyz[i, j, 3, e] = x[3]
            end
        end
    end

    # ============================================================
    # Rebuild ClimaCore local_geometry on Comodo geometry
    # ============================================================

    grid = ClimaCore.Spaces.grid(space)
    global_geometry = getfield(grid, :global_geometry)
    native_local_geometry = getfield(grid, :local_geometry)
    local_geometry_comodo = copy(native_local_geometry)

    FT = Float64
    CoordType2D = ClimaCore.Geometry.LatLongPoint{FT}
    AIdx = ClimaCore.Geometry.coordinate_axis(CoordType2D)

    function comodo_dxyz_dξ(comodo_xyz, elem, i, j, quad_points, Dmat)
        d1 = zeros(FT, 3)
        for ip in eachindex(quad_points)
            x = comodo_xyz_at_ξη_from_face(Vn, Fn, face_match[elem], quad_points[ip], quad_points[j])
            d1 .+= Dmat[i, ip] .* x
        end

        d2 = zeros(FT, 3)
        for jp in eachindex(quad_points)
            x = comodo_xyz_at_ξη_from_face(Vn, Fn, face_match[elem], quad_points[i], quad_points[jp])
            d2 .+= Dmat[j, jp] .* x
        end

        return hcat(d1, d2)
    end

    for elem in 1:nelems
        local_geometry_slab = ClimaCore.DataLayouts.slab(local_geometry_comodo, elem)

        for i in 1:Nq, j in 1:Nq
            x = vec(comodo_xyz[i, j, :, elem])

            xcart = ClimaCore.Geometry.Cartesian123Point(x[1], x[2], x[3])
            u = ClimaCore.Geometry.LatLongPoint(xcart, global_geometry)

            dxyz = comodo_dxyz_dξ(comodo_xyz, elem, i, j, quad_points, Dmat)

            ∂x∂ξ = ClimaCore.Geometry.AxisTensor(
                (
                    ClimaCore.Geometry.Cartesian123Axis(),
                    ClimaCore.Geometry.CovariantAxis{AIdx}(),
                ),
                dxyz,
            )

            G = ClimaCore.Geometry.local_to_cartesian(global_geometry, u)
            ∂u∂ξ = ClimaCore.Geometry.project(
                ClimaCore.Geometry.LocalAxis{AIdx}(),
                G' * ∂x∂ξ,
            )

            J = det(ClimaCore.Geometry.components(∂u∂ξ))
            WJ = J * quad_weights[i] * quad_weights[j]

            local_geometry_slab[ClimaCore.DataLayouts.slab_index(i, j)] =
                ClimaCore.Geometry.LocalGeometry(u, J, WJ, ∂u∂ξ)
        end
    end

    geom_comodo_exact = ClimaCore.Fields.Field(local_geometry_comodo, space)

    # ============================================================
    # 3D mesh visualization (clean)
    # ============================================================

    fig_mesh = Figure(size = (1200, 600))

    ax1 = Axis3(fig_mesh[1, 1], aspect = :data, title = "ClimaCore Mesh ($(mapping_type))")
    ax2 = Axis3(fig_mesh[1, 2], aspect = :data, title = "Comodo Mesh")

    # --- ClimaCore mesh ---
    for e in 1:nelems
        Xc = clima_xyz[:, :, 1, e]
        Yc = clima_xyz[:, :, 2, e]
        Zc = clima_xyz[:, :, 3, e]

        for j in 1:Nq
            lines!(ax1, Xc[:, j], Yc[:, j], Zc[:, j], color = (:gray, 0.5), linewidth = 1)
        end
        for i in 1:Nq
            lines!(ax1, Xc[i, :], Yc[i, :], Zc[i, :], color = (:gray, 0.5), linewidth = 1)
        end

        scatter!(ax1, vec(Xc), vec(Yc), vec(Zc), markersize = 5, color = :blue)
    end

    # --- Comodo mesh ---
    vx_comodo = [p[1] for p in Vn]
    vy_comodo = [p[2] for p in Vn]
    vz_comodo = [p[3] for p in Vn]
    for f in eachindex(Fn)
        i1, i2, i3, i4 = Tuple(Fn[f])
        cyc = [i1, i2, i3, i4, i1]
        lines!(
            ax2,
            [vx[k] for k in cyc],
            [vy[k] for k in cyc],
            [vz[k] for k in cyc],
            color = (:gray, 0.5),
            linewidth = 1
        )
    end

    scatter!(ax2, vx_comodo, vy_comodo, vz_comodo, markersize = 5, color = :blue)

    save("mesh_comparison_$(mapping_type).png", fig_mesh)
    display(fig_mesh)

    # ============================================================
    # J plots + statistics
    # ============================================================

    J_clima  = vec(Array(parent(geom.J)))
    J_comodo = vec(Array(parent(geom_comodo_exact.J)))

    stats_clima  = stats_summary(J_clima)
    stats_comodo = stats_summary(J_comodo)

    println("\n=== J STATISTICS ($(mapping_type)) ===")
    print_stats("ClimaCore ($(mapping_type)):", stats_clima)
    print_stats("Comodo:", stats_comodo)

    figJ = Figure(size = (1200, 400))

    axJ1 = Axis(figJ[1, 1],
        title = "ClimaCore ($(mapping_type)): J vs Point Index",
        xlabel = "Point Index",
        ylabel = "|J determinant|"
    )
    lines!(axJ1, 1:length(J_clima), J_clima, linewidth = 1.5)

    axJ2 = Axis(figJ[1, 2],
        title = "Comodo: J vs Point Index",
        xlabel = "Point Index",
        ylabel = "|J determinant|"
    )
    lines!(axJ2, 1:length(J_comodo), J_comodo, linewidth = 1.5)

    save("J_comparison_$(mapping_type).png", figJ)
    display(figJ)
end