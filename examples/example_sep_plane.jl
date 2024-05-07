using PolyPlanning

N_polys = 3
ego_polys = PolyPlanning.gen_ego_rect()
sep_prob = PolyPlanning.setup_sep_planes(ego_polys;
    T=20,
    dt=.1,
    sides_per_poly=4
);
#polys = PolyPlanning.gen_polys(N_polys)
(P1, P2, P3) = polys

# to display:
PolyPlanning.plot_polys(polys)

x0 = [-3., 2.0, .5, 0, 0, 0];
sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0, P1, P2, P3);

# compare with our method:
#A, b = PolyPlanning.poly_from(x0, sep_prob.angles, sep_prob.lengths)
#self_poly = PolyPlanning.ConvexPolygon2D(A, b)
#ego_polys = [self_poly]

our_prob = PolyPlanning.setup_quick(ego_polys;
    N_polys,
    T=20,
    dt=.1,
    sides_per_poly=4
);
our_sol = PolyPlanning.solve_quick(our_prob, x0, polys);