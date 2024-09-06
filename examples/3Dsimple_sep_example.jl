using PolyPlanning

r=rand
mrp = (ones(3)*0.8+r(3)*0.4)/sqrt(3) # tan(θ/4)∈[0.8, 1.2] θ∈[2.70, 3.50]
# e, θ = PolyPlanning.axis_angle_from_mrp(mrp)
# err = mrp - PolyPlanning.mrp_from_axis_angle(e, θ)
# if norm(err)>1e-4
#     @warn err
#     println(e)
#     println(rad2deg(θ))
# end
trans =zeros(3) + [5,2,3] + r(3)
x0 = [trans; mrp; zeros(6)]

Ve = [[-r(), -r(), -r()], [r(), -r(), -r()], [0, r(), -r()], [0, 0, 5r()]]
Vo = [[-1-r(), -1-r(), -1-r()], [1+r(), -1-r(), -1-r()], [0, 1+r(), -1-r()], [0, 0, 5r()]]

x0 = [ 5.129113794516678,2.515060674726819,3.2043899776388827,0.38173923415793265,0.7635567015039818,0.5014663742804137,0.0,0.0,0.0,0.0,0.0,0.0]
Ve = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]
Vo = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]

Pe = PolyPlanning.ConvexPolygon3D(Ve)
Po = PolyPlanning.ConvexPolygon3D(Vo)
ego_polys = [Pe]
obs_polys = [Po]
R_cost = 1e-3 * PolyPlanning.I(6)
R_cost[4:6, 4:6] = R_cost[4:6, 4:6] / 100.0

sep_prob = PolyPlanning.setup_sep_planes_3d(
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
)

sep_sol = PolyPlanning.solve_prob_sep_planes_3d(sep_prob, x0; is_displaying=true, sleep_duration=0.1)
