using PolyPlanning
ego_L = PolyPlanning.gen_ego_L();
prob = PolyPlanning.setup_quick(
    ego_L;
    T=50,
    dt=0.2,
    Q=0.01 * [1.0 0; 0 1],
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
    N_polys=2
);
gap_polys = PolyPlanning.gen_gap();
x0 = [-5.5, 0, pi / 2, 0, 0, 0];
sol = PolyPlanning.solve_quick(prob, x0, gap_polys);