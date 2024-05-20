using PolyPlanning

N_polys = 3
polys = PolyPlanning.gen_polys(N_polys);
# to display:
PolyPlanning.plot_polys(polys);
x0 = [.0, 1.0, 0.1, 0, 0, 0];

ego_rect = PolyPlanning.gen_ego_rect()
(P1, P2, P3) = polys;

sep_prob = PolyPlanning.setup_sep_planes(
    ego_rect,
    polys;
    T=20,
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
    N_polys=3
);

sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0);

# compare with our method:
#our_prob = PolyPlanning.setup_quick(
#    ego_rect;
#    T=20,
#    dt=0.1,
#    L=1.0,
#    Q=0.0 * [1.0 0; 0 1],
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
#    N_polys
#);
#our_sol = PolyPlanning.solve_quick(our_prob, x0, polys);