using PolyPlanning
using Polyhedra
using LinearAlgebra
using SparseArrays
using GLMakie
using Symbolics

mrp = ([4,2,1]+r(3))/4
# e, θ = PolyPlanning.axis_angle_from_mrp(mrp)
# err = mrp - PolyPlanning.mrp_from_axis_angle(e, θ)
# if norm(err)>1e-4
#     @warn err
#     println(e)
#     println(rad2deg(θ))
# end
trans =zeros(3) + [8,0,0] + r(3)
x0 = [trans; mrp; zeros(6)]

r=rand
Ve1 = [[-1.0,-1,-1], [1,-1,-1], [-1,1,-1], [1,1,-1], [-1,-1,0.8-1], [1,-1,0.8-1], [-1,1,0.8-1], [1,1,0.8-1]]
Ve2 = [[1.0,-1,-1], [1.8,-1,-1], [1,1,-1], [1.8,1,-1], [1.0,-1,2-1], [1.8,-1,2-1], [1,1,2-1], [1.8,1,2-1]]
Vo1 = [[1.0,-5,1],[1,5,1],[1,-5,5],[1,5,5], [2.0,-5,1],[2,5,1],[2,-5,5],[2,5,5]]
Vo2 = [[1.0,-5,-1],[1,5,-1],[1,-5,-5],[1,5,-5], [2.0,-5,-1],[2,5,-1],[2,-5,-5],[2,5,-5]]

Pe1 = PolyPlanning.ConvexPolygon3D(Ve1)
Pe2 = PolyPlanning.ConvexPolygon3D(Ve2)
Po1 = PolyPlanning.ConvexPolygon3D(Vo1)
Po2 = PolyPlanning.ConvexPolygon3D(Vo2)
ego_polys = [Pe1, Pe2]
obs_polys = [Po1, Po2]

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

our_sol = PolyPlanning.solve_nonsmooth_3d(nonsmooth_prob, x0; is_displaying=true, sleep_duration=0.1)