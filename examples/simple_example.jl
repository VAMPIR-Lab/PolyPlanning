using PolyPlanning

#x0 = [5.0, 0.0, 0.1, 0, 0, 0];
x0 = [5.0, 0.0, -1.1, 0, 0, 0]; #mcp no progress
# x0 = [5.0, 0.0, -1.0, 0, 0, 0]; #mcp is solved but seems not convergent
# x0 = [5.0, 0.0, -.5, 0, 0, 0]; #mcp is solved but gets stuck when one corner of ego in contact with the obstacle
# x0 = [5.0, 0.0, -.3, 0, 0, 0]; #mcp is solved but gets stuck when one corner of ego in contact with the obstacle
# x0 = [5.0, 0.0, -.0, 0, 0, 0]; #mcp is solved but gets stuck when short edge of ego in contact with the obstacle
# x0 = [5.0, 0.0, .1, 0, 0, 0]; #mcp is solved and reaches the goal
# x0 = [5.0, 0.0, .2, 0, 0, 0]; #mcp is solved and reaches the goal
# x0 = [5.0, 0.0, .3, 0, 0, 0]; #mcp is solved but gets stuck when one corner of ego in contact with the obstacle
obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);
ego_rect = PolyPlanning.gen_ego_rect(; a=0.5, b=1.0);
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;
is_newsd=true

our_prob = PolyPlanning.setup_quick(
    ego_rect;
    T=20,
    dt=0.2,
    Q=0.0 * PolyPlanning.I(2), # disabled final cost
    q=[0, 0.0],
    Rf,
    Qf=5e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
    n_obs=length(obs_polys),
    is_newsd=is_newsd
);

our_sol = PolyPlanning.solve_quick(our_prob, x0, obs_polys; is_displaying=true, sleep_duration = 1.5, is_newsd=is_newsd)


