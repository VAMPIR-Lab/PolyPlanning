using PolyPlanning

x0 = [5.0, 0.0, π / 2 + .1, 0, 0, 0];
obs_polys = PolyPlanning.gen_rect_obs(; a=0.5, b=2.0);
ego_polys = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);

R_cost = 1e-3 * PolyPlanning.I(3);
R_cost[3, 3] = R_cost[3, 3] / 100.0;

prob = PolyPlanning.setup_nonsmooth(
    ego_polys,
    obs_polys;
    T=20,
    dt=.2,
    R_cost,
    Q_cost=2e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π
)

PolyPlanning.solve_nonsmooth(prob, x0; is_displaying=true, sleep_duration=0.25)