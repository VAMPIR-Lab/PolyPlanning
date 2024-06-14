using PolyPlanning

#x0 = [5.0, 0.0, rand(), 0, 0, 0];
#x0 = [5.0, 0.0, -π / 2 + π * rand(), 0, 0, 0];
obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);
#obs_polys = PolyPlanning.gen_simple_obs();
#ego_rect = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);
#Rf = 1e-3 * PolyPlanning.I(3);
ego_rect = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;

our_prob = PolyPlanning.setup_quick(
    ego_rect;
    T=20,
    dt=0.2,
    Rf,
    Qf=2e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
    n_obs=length(obs_polys)
);

#x0 = [5.0, 0.0, -π / 2 + π * rand(), 0, 0, 0];
#x0 = [5.0, 0.0, 0.006216649582387657, 0, 0, 0];
x0 = [5.0, 0.0, -1.2633343742220948, 0, 0, 0];
#x0 = [5.0, 0.0, π/2+, 0, 0, 0];
our_sol = PolyPlanning.solve_quick(our_prob, x0, obs_polys; is_displaying=true, sleep_duration=0.25)


