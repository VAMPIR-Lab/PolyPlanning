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
trans =zeros(3) + [5,2,3] + r(3)
x0 = [trans; mrp; zeros(6)]

r=rand
Ve = [[-r(), -r(), -r()], [r(), -r(), -r()], [0, r(), -r()], [0, 0, 5r()]]
Vo1 = [[-1-r(), -1-r(), -1-r()], [1+r(), -1-r(), -1-r()], [0, 1+r(), -1-r()], [0, 0, 5r()]]
Vo2 = [[-1-r(), -1-r(), 1+r()], [1+r(), -1-r(), 1+r()], [0, 1+r(), 1+r()], [0, 0, -5r()]]
Vo3 = [[4r(),4r(),4r()], [4r(),4r(),4r()], [4r(),4r(),4r()], [4r(),4r(),4r()]]


# success example
x0=[5.650530044055372,2.543840168058619,3.07638181761498,1.026193123994465,0.6354979144725007,0.4713567932520325,0.0,0.0,0.0,0.0,0.0,0.0]
Ve = [[-0.31822548564254016, -0.3581025253582616, -0.17516160728486496],
[0.5329421287163413, -0.6858966368162803, -0.49184278399090586],
[0.0, 0.46256837919661886, -0.3630302271283352], [0.0, 0.0, 0.8325626114258877]]
Vo1=[[-1.4124579489604305, -1.7408834343392958, -1.183853525550069],
[1.5782362174819262, -1.34633918610923, -1.1360703212170882],
[0.0, 1.8319664436476089, -1.4708332611212196] ,
[0.0, 0.0, 4.520406526968425]
]
Vo2=[[-1.8090249024122276, -1.837362900170374, 1.4926843679435198],
[1.6210436341615457, 1.6172860789861354, 1.7076346960493693],
[0.0, 1.118147950193499, 1.462411504373387]    ,
[0.0, 0.0, -1.9556185996889925]
]
Vo3=[ [3.3718050377688256, 2.533919251941978, 1.062872982009209],
[1.736542856536936, 2.5761102539537135, 2.277017652819808],
[0.7677146554669978, 1.5569621989157247, 1.1063109456932385],
[3.1306454412167373, 3.9037170982605693, 3.3219174189108562]
]

Pe = PolyPlanning.ConvexPolygon3D(Ve)
Po1 = PolyPlanning.ConvexPolygon3D(Vo1)
Po2 = PolyPlanning.ConvexPolygon3D(Vo2)
Po3 = PolyPlanning.ConvexPolygon3D(Vo3)
ego_polys = [Pe]
obs_polys = [Po1, Po2, Po3]


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

our_sol = PolyPlanning.solve_nonsmooth_3d(nonsmooth_prob, x0; is_displaying=true, sleep_duration=0.01)