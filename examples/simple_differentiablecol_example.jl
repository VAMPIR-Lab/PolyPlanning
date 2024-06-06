using PolyPlanning

x0 = [5.0, 0.0, 0.1, 0, 0, 0];
obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);
ego_rect = PolyPlanning.gen_ego_rect(; a=0.5, b=1.0);
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;

our_prob = PolyPlanning.setup_differentiablecol(
    ego_rect;
    T=1,
    dt=0.2,
    Rf,
    Qf=5e-3 * PolyPlanning.I(2),
    u1_max=10*10.0,
    u2_max=10*10.0,
    u3_max=10*π,
    n_obs=length(obs_polys)
);

our_sol = PolyPlanning.solve_differentiablecol(our_prob, x0, obs_polys; is_displaying=true)
