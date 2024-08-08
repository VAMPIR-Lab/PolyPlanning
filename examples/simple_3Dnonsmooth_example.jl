using PolyPlanning
using Polyhedra
using LinearAlgebra
using SparseArrays
using GLMakie
using Symbolics

# U shape
# V1 = [[0.0, 0, 0], [0, 0, 1], [1, 0, 0], [1, 0, 1], [0, 1, 0], [0, 1, 1], [1, 1, 0], [1, 1, 1]]
# V2 = [[2, 0.0, 0], [2, 0, 1], [3, 0, 0], [3, 0, 1], [2, 1, 0], [2, 1, 1], [3, 1, 0], [3, 1, 1]]
# V3 = [[0, 1, 0.0], [0, 1, 1], [3, 1, 0], [3, 1, 1], [0, 2, 0], [0, 2, 1], [3, 2, 0], [3, 2, 1]]
# P1 = PolyPlanning.ConvexPolygon3D(V1)
# P2 = PolyPlanning.ConvexPolygon3D(V2)
# P3 = PolyPlanning.ConvexPolygon3D(V3)
# ego_polys = [P1, P2, P3]
# obs_polys = [P1, P2, P3]

# Ve = [[.25, -2, -1], [.25, 2, -1], [-.25, 2, -1], [-.25, -2, -1], [-.5, -.5, 1], [.5, -.5, 1], [-.5, .5, 1], [.5, .5, 1]]
# Vo = [[-1.0, -1, -1], [1, -1, -1], [-1, 1, -1], [1, 1, -1], [0, 0, 5]]
# Ve2 = -1*Ve
# Pe2 = PolyPlanning.ConvexPolygon3D(Ve2)
# Vo2 = -1*Vo
# Po2 = PolyPlanning.ConvexPolygon3D(Vo2)

Ve = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]
Vo = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]
Pe = PolyPlanning.ConvexPolygon3D(Ve)
Po = PolyPlanning.ConvexPolygon3D(Vo)
ego_polys = [Pe]
obs_polys = [Po]

mrp = [4,2,1]/4
e, θ = PolyPlanning.axis_angle_from_mrp(mrp)
err = mrp - PolyPlanning.mrp_from_axis_angle(e, θ)
# if norm(err)>1e-4
#     @warn err
#     println(e)
#     println(rad2deg(θ))
# end
trans =zeros(3)+ [5,2,3]
x0 = [trans; mrp; zeros(6)]

# Ve = [[1, 1], [2, 1], [2.5, 1.4], [2.5, 2], [1.5, 1.8], [.8, 1.2]]
# Ve2 = [[2, 1], [3, -2], [2., -2], [3., 0], [.8, -1.2]]
# ego_polys = [PolyPlanning.ConvexPolygon2D(Ve), PolyPlanning.ConvexPolygon2D(Ve2)]
# Vo = [[.25, -2], [.25, 2], [-.25, 2], [-.25, -2]]
# Vo2 = [[0, 0.0], [1, -1], [1, 1]]
# obs_polys = [PolyPlanning.ConvexPolygon2D(Vo), PolyPlanning.ConvexPolygon2D(Vo2)]

R_cost = 1e-3 * PolyPlanning.I(6)
R_cost[4:6, 4:6] = R_cost[4:6, 4:6] / 100.0
nonsmooth_prob = PolyPlanning.setup_nonsmooth_3d(
    ego_polys,
    obs_polys;
    T=20,
    dt=.2,
    R_cost,
    Q_cost=2e-3 * PolyPlanning.I(3),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=10.0,
    u4_max=π,
    u5_max=π,
    u6_max=π,
	n_sd_slots=15
)

our_sol = PolyPlanning.solve_nonsmooth_3d(nonsmooth_prob, x0; is_displaying=true)#, sleep_duration=0.01)