using PolyPlanning

x0 = [5.0, 0.0, 0.1, 0, 0, 0];
obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);
ego_rect = PolyPlanning.gen_ego_rect(; a=0.5, b=1.0);
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;

kkt_prob = PolyPlanning.setup_direct_kkt(
    ego_rect,
    obs_polys;
    T=20,
    dt=0.2,
    Rf,
    Qf=5e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π
);

kkt_sol = PolyPlanning.solve_prob_direct_kkt(kkt_prob, x0; is_displaying=true)
