# metrics: 
# flat wall calibration
# random initial conditions
# random walls
# random ego shape? (later)
# time to solve
# success from mcp status
# success (from constraints)
# x dist

using PolyPlanning

n_walls = 10
n_x0s = 10
n_sides = 4
n_obs = 4
T = 20
dt = 0.2
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 10.0;
Qf = 1e-2 * PolyPlanning.I(2),
p1_max = 50.0
p2_min = -50.0
u1_max = 10.0
u2_max = 10.0
u3_max = π

# calibration
#obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);
#ego_rect = PolyPlanning.gen_ego_rect(; l_multip=2.0);
#x0 = [2.0, 2.0, 0.1, 0, 0, 0]
#prob = PolyPlanning.setup_quick(
#    ego_rect;
#    T,
#    dt,
#    Q=0.0 * [1.0 0; 0 1], # disabled final cost
#    q=[0, 0.0],
#    Rf,
#    Qf,
#    p1_max,
#    p2_min,
#    u1_max,
#    u2_max,
#    u3_max,
#    sides_per_poly=4,
#    derivs_per_sd=4,
#    derivs_per_fv=4,
#    N_polys=length(obs_polys)
#);
#sol = PolyPlanning.solve_quick(prob, x0, obs_polys)

# example
obs_polys = PolyPlanning.gen_packing_wall(n_obs, n_sides; w=1.5, l=2.0, max_overlap=0.0);
PolyPlanning.plot_polys(obs_polys);
ego_rect = PolyPlanning.gen_ego_rect(; l_multip=2.0);
x0 = [2.0, 2.0, 0.1, 0, 0, 0]

#PolyPlanning.plot_polys(obs_polys);
#walls = map(1:n_maps) do i
#    PolyPlanning.gen_packing_wall(n_obs, n_sides; width=1.0, length=2.0)
#    PolyPlanning.plot_polys(obs_polys)
#end

#x0s = map(1:n_x0s) do i
#    [1.5, 0.0, 0.1, 0, 0, 0]
#end

prob = PolyPlanning.setup_quick(
    ego_rect;
    T,
    dt,
    Q=0.0 * [1.0 0; 0 1], # disabled final cost
    q=[0, 0.0],
    Rf,
    Qf,
    p1_max,
    p2_min,
    u1_max,
    u2_max,
    u3_max,
    sides_per_poly=4,
    derivs_per_sd=4,
    derivs_per_fv=4,
    N_polys=length(obs_polys)
);

sol = PolyPlanning.solve_quick(prob, x0, obs_polys)

#direct_prob = PolyPlanning.setup_direct_kkt(
#    ego_rect,
#    obs_polys;
#    T,
#    dt,
#    R,
#    p1_max,
#    p2_min,
#    u1_max,
#    u2_max,
#    u3_max
#);

#direct_sol = PolyPlanning.solve_prob_direct_kkt(direct_prob, x0);

#sep_prob = PolyPlanning.setup_sep_planes(
#    ego_rect,
#    obs_polys;
#    T,
#    dt,
#    R,
#    p1_max,
#    p2_min,
#    u1_max,
#    u2_max,
#    u3_max
#);

#sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0);