using PolyPlanning

x0 = [2.5, 0.0, 0.1, 0, 0, 0]
#obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);
obs_polys = PolyPlanning.gen_packing_wall(2, 4);
ego_rect = PolyPlanning.gen_ego_rect(; l_multip=2.0);
T = 20
Rf = 1e-3 * PolyPlanning.I(3)
Rf[3, 3] = Rf[3, 3] / 10.0

sep_prob = PolyPlanning.setup_sep_planes(
    ego_rect,
    obs_polys;
    T,
    dt=0.2,
    Rf,
    Qf=1e-2 * PolyPlanning.I(2),
    u1_max = 1.0,
    u2_max=1.0,
    u3_max=π
)

sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0; is_displaying=false)

(fig, update_fig) = PolyPlanning.visualize_sep_planes(x0, T, ego_rect, obs_polys)
update_fig(sep_sol.θ)
display(fig)
