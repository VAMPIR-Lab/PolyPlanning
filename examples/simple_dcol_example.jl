using PolyPlanning

x0 = [5.0, 0.0, π/2 + .1, 0, 0, 0];
#obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);
#ego_rect = PolyPlanning.gen_ego_rect(; a=0.5, b=1.0);
obs_polys = PolyPlanning.gen_rect_obs(; a=0.5, b=2.0);
ego_polys = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;

dcol_prob = PolyPlanning.setup_dcol(
    ego_polys;
    T=20,
    dt=0.2,
    Rf,
    Qf=1e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
    n_obs=length(obs_polys)
);

dcol_sol = PolyPlanning.solve_dcol(dcol_prob, x0, obs_polys; is_displaying=true, sleep_duration=0.1)
