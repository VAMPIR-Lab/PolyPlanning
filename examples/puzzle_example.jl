using PolyPlanning

x0 = [5.0, -3.0, -π/4, 0.0, 0.0, 0.0]; 

#ego_length = 2.0
#ego_width = 0.5

obs_polys = PolyPlanning. gen_T_obstacle(; base_width=5.0, width=1.5, length=1.0)
ego_polys = PolyPlanning.gen_ego_U();

R_cost = 1e-3 * PolyPlanning.I(3);
R_cost[3, 3] = R_cost[3, 3] / 100.0;

nonsmooth_prob = PolyPlanning.setup_nonsmooth(
    ego_polys,
    obs_polys;
    T=20,
    dt=0.2,
    R_cost,
    Q_cost=2e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
    n_sd_slots=8
)

PolyPlanning.solve_nonsmooth(nonsmooth_prob, x0; is_displaying=true, sleep_duration=0.01)