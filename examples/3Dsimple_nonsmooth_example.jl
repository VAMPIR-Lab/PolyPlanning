using PolyPlanning
using Polyhedra
using LinearAlgebra
using SparseArrays
using GLMakie
using Symbolics

r=rand
# mrp = ([4,2,1]+r(3))/4 # seems to get a higher success rate
# mrp = (ones(3)*0.8+r(3)*0.4)/sqrt(3) # tan(θ/4)∈[0.8, 1.2] θ∈[2.70, 3.50]
axis = rand(3)
mrp = axis / norm(axis) * (0.8+0.4*r()) # norm(mrp)∈[0.8, 1.2]
# e, θ = PolyPlanning.axis_angle_from_mrp(mrp)
# err = mrp - PolyPlanning.mrp_from_axis_angle(e, θ)
# if norm(err)>1e-4
#     @warn err
#     println(e)
#     println(rad2deg(θ))
# end
trans =zeros(3) + [5,2,3] + r(3)
x0 = [trans; mrp; zeros(6)]

Ve = [[-1-r(), -1-r(), -1-r()], [1+r(), -1-r(), -1-r()], [0, 1+r(), -1-r()], [0, 0, 1+r()]]
Vo = [[-1-r(), -1-r(), -1-r()], [1+r(), -1-r(), -1-r()], [0, 1+r(), -1-r()], [0, 0, 5r()]]

# success example
# Ve = [[0.3557744694359173, 0.7163815945257631, 0.6655002776158867],
# [-0.3359690949338091, -0.46250510373333853, -0.44019465179445494],
# [1.3876161809838614, 0.06523144992337127, 0.9657548361522916],
# [-0.055131582061873186, -0.9003721133659086, -1.3374867091607676]]
# Vo = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]

# success example
# Ve = [ [-0.8107863027036382, -0.36669554366980484, -0.4491413389259594],
# [0.5228937962319313, -0.6877026417799523, -0.5654583041992148],
# [0.0, 0.9823155313950472, -0.3550657247336866] ,
# [0.0, 0.0, 3.5857541379504765]
# ]
# Vo = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]

# success example
# Ve = [[-1.0, -1, -1.1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]
# Vo = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]

# failure example
# Ve = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]
# Vo = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]


Pe = PolyPlanning.ConvexPolygon3D(Ve)
Po = PolyPlanning.ConvexPolygon3D(Vo)
ego_polys = [Pe]
obs_polys = [Po]


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