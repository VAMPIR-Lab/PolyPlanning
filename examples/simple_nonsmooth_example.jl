using PolyPlanning
ego_rect = PolyPlanning.gen_ego_rect();

n_obs = 1
obs_polys = PolyPlanning.gen_polys(n_obs, side_length=4); PolyPlanning.plot_polys(obs_polys);
x0 = [.5, -1, .1, 0, 0, 0];

prob = PolyPlanning.setup_quick(
    ego_rect;
    T=20,
    dt=0.2,
    Q=0.01 * [1.0 0; 0 1],
    q=[0, 0.0],
    R=0.01 * PolyPlanning.I(3),
    p1_max=500.0,
    p2_min=-500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=1*π / 4,
    sides_per_poly=4,
    derivs_per_sd=4,
    derivs_per_fv=4,
    N_polys=length(obs_polys)
);

sol = PolyPlanning.solve_quick(prob, x0, obs_polys)