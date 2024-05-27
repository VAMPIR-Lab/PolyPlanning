using PolyPlanning

n_obs = 3
obs_polys = PolyPlanning.gen_polys(n_obs, side_length=4); PolyPlanning.plot_polys(obs_polys);
x0 = [-3., 1, .1, 0, 0, 0];

# or you could try
#obs_polys = PolyPlanning.gen_gap()
#ego_rect = PolyPlanning.gen_ego_L()
#x0 = [-5.5, 0, pi / 2, 0, 0, 0];
#ego_rect = PolyPlanning.gen_ego_U()
#x0 = [-5,.2,pi/2+.2,0,0,0]

sep_prob = PolyPlanning.setup_sep_planes(
    ego_rect,
    obs_polys;
    T=20,
    dt=0.2,
    R=0.01 * PolyPlanning.I(3),
    p1_max=500.0,
    p2_min=-500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=1*π / 4
);

sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0);

# compare with our method:
#our_prob = PolyPlanning.setup_quick(
#    ego_rect;
#    T=20,
#    dt=0.2,
#    Q= 0*[1.0 0; 0 1],
#    q=[0, 0.0],
#    R=0.01 * PolyPlanning.I(3),
#    p1_max=500.0,
#    p2_min=-500.0,
#    u1_max=1.0,
#    u2_max=1.0,
#    u3_max=π / 4,
#    sides_per_poly=4,
#    derivs_per_sd=4,
#    derivs_per_fv=4,
#    N_polys=length(obs_polys)
#);
#our_sol = PolyPlanning.solve_quick(our_prob, x0, obs_polys);