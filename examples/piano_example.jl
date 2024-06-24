using PolyPlanning

#x0 = [3.2255779100589113, -3.769459974984243, 0.0, 0.0, 0.0, 0.0]; # success
x0 = [3.493367162884191, -3.830529204368667, 0.0, 0.0, 0.0, 0.0]; # fail
pre_L_length_base = 5.0
post_L_length = 3.0
ego_length = 2.0
ego_width = 0.5
width = sqrt((ego_length / 2)^2 + ego_width^2) + 0.2
pre_L_length = pre_L_length_base - width / 2
obs_polys = PolyPlanning.gen_L_corridor(; width, pre_L_length, post_L_length)
ego_polys = PolyPlanning.gen_ego_rect(; a=ego_width, b=ego_length);

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
    n_sd_slots=4
)

PolyPlanning.solve_nonsmooth(nonsmooth_prob, x0; is_displaying=true, sleep_duration=0.01)