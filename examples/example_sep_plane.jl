using PolyPlanning

N_polys = 3
#polys = PolyPlanning.gen_polys(N_polys);
# to display:
PolyPlanning.plot_polys(polys);
x0 = [0.0, -5.0, 0.1, 0, 0, 0];

ego_polys = PolyPlanning.gen_ego_rect()
sep_prob = PolyPlanning.setup_sep_planes(ego_polys;
    T=5,
    dt=0.5,
    sides_per_poly=4
);
(P1, P2, P3) = polys;
sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0, P1, P2, P3);

# compare with our method:
our_prob = PolyPlanning.setup_quick(ego_polys;
    N_polys,
    T=5,
    dt=0.5,
    sides_per_poly=4
);
our_sol = PolyPlanning.solve_quick(our_prob, x0, polys);