using PolyPlanning

x0 = [5.0, 0.0, π / 2 + .1, 0, 0, 0];
obs_polys = PolyPlanning.gen_rect_obs(; a=0.5, b=2.0);
ego_polys = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);

# Ve = [[1, 1], [2, 1], [2.5, 1.4], [2.5, 2], [1.5, 1.8], [.8, 1.2]]
# Ve2 = [[2, 1], [3, -2], [2., -2], [3., 0], [.8, -1.2]]
# ego_polys = [PolyPlanning.ConvexPolygon2D(Ve), PolyPlanning.ConvexPolygon2D(Ve2)]
# Vo = [[.25, -2], [.25, 2], [-.25, 2], [-.25, -2]]
# Vo2 = [[0, 0.0], [1, -1], [1, 1]]
# obs_polys = [PolyPlanning.ConvexPolygon2D(Vo), PolyPlanning.ConvexPolygon2D(Vo2)]

R_cost = 1e-3 * PolyPlanning.I(3);
R_cost[3, 3] = R_cost[3, 3] / 100.0;

nonsmooth_prob = PolyPlanning.setup_nonsmooth(
    ego_polys,
    obs_polys;
    T=20,
    dt=.2,
    R_cost,
    Q_cost=2e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
	n_sd_slots=4
)

PolyPlanning.solve_nonsmooth(nonsmooth_prob, x0; is_displaying=true, sleep_duration=0.25)