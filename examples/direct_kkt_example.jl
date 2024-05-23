using PolyPlanning

n_obs = 1
#obs_polys = PolyPlanning.gen_polys(n_obs, side_length=4); PolyPlanning.plot_polys(obs_polys);
x0 = [-1.0, 4.0, 0.1, 0, 0, 0];

# need support for multiple ego polys
ego_rect = PolyPlanning.gen_ego_rect()

direct_prob = PolyPlanning.setup_direct_kkt(
    ego_rect,
    obs_polys;
    T=10,
    dt=0.1,
    R=0.01 * PolyPlanning.I(3),
    p1_max=500.0,
    p2_min=-500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=1*π / 4,
    sides_per_obs=length(obs_polys[1].b),
    sides_per_ego=length(ego_rect[1].b)
);

direct_sol = PolyPlanning.solve_prob_direct_kkt(direct_prob, x0);

# compare with our method:
our_prob = PolyPlanning.setup_quick(
    ego_rect;
    T=10,
    dt=0.1,
    L=1.0,
    Q=0.0 * [1.0 0; 0 1],
    q=[0, 0.0],
    R=0.01 * PolyPlanning.I(3),
    p1_max=500.0,
    p2_min=-500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=π / 4,
    sides_per_poly=4,
    derivs_per_sd=4,
    derivs_per_fv=4,
    N_polys
);
our_sol = PolyPlanning.solve_quick(our_prob, x0, polys);