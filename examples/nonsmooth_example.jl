using PolyPlanning

x0 = [6.0, 2.0, 0.1, 0, 0, 0]
obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);;
ego_rect = PolyPlanning.gen_ego_rect(; l_multip=2.0);
T = 20
Rf = 1e-3 * PolyPlanning.I(3)
Rf[3, 3] = Rf[3, 3] / 10.0

our_prob = PolyPlanning.setup_quick(
    ego_rect;
    T,
    dt=0.2,
    Q=0.0 * [1.0 0; 0 1], # disabled final cost
    q=[0, 0.0],
    Rf,
    Qf=1e-2 * PolyPlanning.I(2),
    u1_max=1.0,
    u2_max=1.0,
    u3_max=π,
    n_obs=length(obs_polys)
);

our_sol = PolyPlanning.solve_quick(our_prob, x0, obs_polys; is_displaying=false)

(fig, update_fig) = PolyPlanning.visualize_quick(x0, T, ego_rect, obs_polys)
update_fig(our_sol.θ)
display(fig)

