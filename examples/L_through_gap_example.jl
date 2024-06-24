using PolyPlanning

x0 = [-5.5, 0, pi / 2, 0, 0, 0];
ego_polys = PolyPlanning.gen_ego_L();
gap_polys = PolyPlanning.gen_gap();
R_cost = 1e-3 * PolyPlanning.I(3);
R_cost[3, 3] = R_cost[3, 3] / 100.0;

nonsmooth_prob = PolyPlanning.setup_nonsmooth(
    ego_polys,
    gap_polys;
    T=20,
    dt=0.2,
    R_cost,
    Q_cost=2e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
    n_sd_slots=2
);
    
our_sol = PolyPlanning.solve_nonsmooth(nonsmooth_prob, x0; is_displaying=true, sleep_duration=0.3)
