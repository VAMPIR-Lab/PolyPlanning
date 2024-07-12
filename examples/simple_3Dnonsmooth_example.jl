using PolyPlanning
using Polyhedra
using LinearAlgebra
using SparseArrays
using GLMakie
using Symbolics


Ve = [[.25, -2, -1], [.25, 2, -1], [-.25, 2, -1], [-.25, -2, -1], [-.5, -.5, 1], [.5, -.5, 1], [-.5, .5, 1], [.5, .5, 1]]
Pe = PolyPlanning.ConvexPolygon3D(Ve)
ego_polys = [Pe]
Vo = [[-1.0, -1, -1], [1, -1, -1], [-1, 1, -1], [1, 1, -1], [0, 0, 5]]
Po = PolyPlanning.ConvexPolygon3D(Vo)
obs_polys = [Po]

mrp = [1,2,3]
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

R_cost = 1e-3 * PolyPlanning.I(3);
R_cost[3, 3] = R_cost[3, 3] / 100.0;

nonsmooth_prob = PolyPlanning.setup_nonsmooth(
    ego_polys,
    obs_polys;
    T=20,
    dt=.2,
    R_cost,
    Q_cost=2e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
	n_sd_slots=4
)

PolyPlanning.solve_nonsmooth(nonsmooth_prob, x0; is_displaying=true, sleep_duration=0.25)