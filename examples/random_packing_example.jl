using PolyPlanning

x0 = [5.0, 0.0, 0.1, 0, 0, 0];
#x0 = x0s[8] # too long
#x0 = x0s[1] # fail
obs_polys = PolyPlanning.gen_packing_wall(4, 4; w=5.0, l=5.0, max_overlap=0.0);
#obs_polys = maps[1]
#obs_polys = [obs_polys[4],obs_polys[2],obs_polys[3]]
#obs_polys = [obs_polys[2],obs_polys[3]]
#PolyPlanning.plot_polys(obs_polys);
ego_rect = PolyPlanning.gen_ego_rect(; a=0.5, b=1.0);
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;

our_prob = PolyPlanning.setup_nonsmooth(
    ego_polys,
    obs_polys;
    T=20,
    dt=0.2,
    R_cost,
    Q_cost=2e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
    n_sd_slots=4
);

our_sol = PolyPlanning.solve_nonsmooth(our_prob, x0; is_displaying=true, sleep_duration=0.01)