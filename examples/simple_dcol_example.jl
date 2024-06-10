using PolyPlanning

x0 = [5.0, 0.0, 0.5, 0, 0, 0];
#obs_polys = PolyPlanning.gen_rect_obs(; a=.25);
obs_polys = PolyPlanning.gen_simple_obs();
ego_rect = PolyPlanning.gen_ego_rect(; a=0.5, b=1.0);
Rf = 0 *1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;

dcol_prob = PolyPlanning.setup_dcol(
    ego_rect;
    T=2,
    dt=0.2,
    Rf,
    Qf=5e-3 * PolyPlanning.I(2),
    u1_max=100 * 10.0,
    u2_max=100 * 10.0,
    u3_max=100 * π,
    n_obs=length(obs_polys)
);

dcol_sol = PolyPlanning.solve_dcol(dcol_prob, x0, obs_polys; is_displaying=true)
