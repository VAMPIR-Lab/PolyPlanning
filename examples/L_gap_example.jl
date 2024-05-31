using PolyPlanning

x0 = [-5.5, 0, pi / 2, 0, 0, 0];
gap_polys = PolyPlanning.gen_gap();
ego_L = PolyPlanning.gen_ego_L();
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;

our_prob = PolyPlanning.setup_quick(
    ego_L;
    T=20,
    dt=0.2,
    Q=1e-3 * PolyPlanning.I(2),
    q=[0, 0.0],
    Rf,
    Qf=5e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
    n_obs=length(gap_polys)
);

our_sol = PolyPlanning.solve_quick(our_prob, x0, gap_polys; is_displaying=true)
