using PolyPlanning

N_polys = 3;
sep_prob = PolyPlanning.setup_sep_planes(;
    T=40,
    dt=0.2,
    sides_per_poly=4
);
polys = PolyPlanning.gen_polys(N_polys);
(P1, P2, P3) = polys;

# to display:
#fig = Figure();
#ax = Axis(fig[1, 1], aspect=DataAspect());
#PolyPlanning.plot!(ax, polys);

x0 = [-2.5, 0, 0, 0, 0, 0];
sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0, P1, P2, P3);

# compare with our method:
A, b = PolyPlanning.poly_from(x0, sep_prob.angles, sep_prob.lengths);
self_poly = PolyPlanning.ConvexPolygon2D(A, b);
ego_polys = [self_poly];

our_prob = PolyPlanning.setup_quick(ego_polys;
    N_polys,
    T=50,
    dt=0.2,
    sides_per_poly=4
);
our_sol = PolyPlanning.solve_quick(our_prob, x0, polys);