# random initial conditions
# random walls
# random ego shape? (later)

# metrics: 
# time to solve
# success from mcp status
# success (from constraints) (later)
# x dist

using PolyPlanning

n_walls = 10
n_x0s = 10
n_sides = 4
n_obs = 4
T = 20
dt = 0.2
Rf = 1e-3 * PolyPlanning.I(3)
Rf[3, 3] = Rf[3, 3] / 10.0
Qf = 1e-2 * PolyPlanning.I(2)
u1_max = 10.0
u2_max = 10.0
u3_max = π
init_x = 2.5
init_y_max = 2.0
n_xu = 9

data_dir = "data"
date_now = PolyPlanning.Dates.format(PolyPlanning.Dates.now(), "YYYY-mm-dd_HHMM")

# saving these 
params = (; n_walls,
    n_x0s,
    n_sides,
    n_obs,
    T,
    dt,
    Rf,
    Qf,
    u1_max,
    u2_max,
    u3_max,
    init_x,
    init_y_max
)
PolyPlanning.jldsave("$(data_dir)/experiment_$(date_now).jld2"; ego_rect, walls, x0s, params)

# example
#obs_polys = PolyPlanning.gen_packing_wall(n_obs, n_sides; w=1.5, l=2.0, max_overlap=0.0);
#PolyPlanning.plot_polys(obs_polys);
#ego_rect = PolyPlanning.gen_ego_rect(; l_multip=2.0);
#x0 = [2.0, 2.0, 0.1, 0, 0, 0]

# gen experiment data

ego_rect = PolyPlanning.gen_ego_rect(; l_multip=2.0);

walls = map(1:n_walls) do i
    PolyPlanning.gen_packing_wall(n_obs, n_sides; w=1.5, l=2.0, max_overlap=0.0)
end

x0s = map(1:n_x0s) do i
    [init_x, -init_y_max + 2 * init_y_max * rand(), -π + 2 * π * rand(), 0, 0, 0]
end

# ours
prob = PolyPlanning.setup_quick(
    ego_rect;
    T,
    dt,
    Q=0.0 * [1.0 0; 0 1], # disabled final cost
    q=[0, 0.0],
    Rf,
    Qf,
    u1_max,
    u2_max,
    u3_max,
    sides_per_poly=4,
    derivs_per_sd=4,
    derivs_per_fv=4,
    n_obs
);

@info "computing our sols"
our_sols = Dict()
start = time()
for (i, wall) in enumerate(walls)
    for (j, x0) in enumerate(x0s)
        res = PolyPlanning.solve_quick(prob, x0, wall; is_displaying=false)
        mcp_success = res.status == PolyPlanning.PATHSolver.MCP_Solved
        time = res.info.total_time
        x_dist = res.z[(T-1)*n_xu+1]

        our_sols[i, j] = (; mcp_success, time, x_dist, res)
    end
end
elapsed = time() - start
@info "our sols elapsed $elapsed"
PolyPlanning.jldsave("$(data_dir)/packing_our_sols_$(date_now).jld2"; our_sols)

# simple processing
#our_idxs = []
#our_times = []
#our_x_dists = []

#for (idx, sol) in our_sols
#    wall_idx = idx[1]
#    x0_idx = idx[2]
#    #Main.@infiltrate
#    if sol.mcp_success
#        push!(our_idxs, idx)
#        push!(our_times, sol.time)
#        push!(our_x_dists, sol.x_dist)
#    end
#end

#i_wall = 2
#j_x0 = 5
#(fig, update_fig) = PolyPlanning.visualize_quick(x0s[j_x0], T, ego_rect, walls[i_wall])
#update_fig(our_sols[i_wall, j_x0].res.θ)
#display(fig)

@info "computing sep plane sols"
sep_plane_sols = Dict()
start = time()
for (i, wall) in enumerate(walls)
    sep_prob = PolyPlanning.setup_sep_planes(
        ego_rect,
        wall;
        T,
        dt,
        Rf,
        Qf,
        u1_max,
        u2_max,
        u3_max
    )
    for (j, x0) in enumerate(x0s)
        res = PolyPlanning.solve_prob_sep_planes(sep_prob, x0; is_displaying=false)
        mcp_success = res.status == PolyPlanning.PATHSolver.MCP_Solved
        time = res.info.total_time
        x_dist = res.z[(T-1)*n_xu+1]

        sep_plane_sols[i, j] = (; mcp_success, time, x_dist, res)
    end
end
elapsed = time() - start
@info "sep plane sols elapsed $elapsed"
PolyPlanning.jldsave("$(data_dir)/packing_sep_plane_sols_$(date_now).jld2"; sep_plane_sols)


#sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0);





@info "computing direct sols"

#sol = PolyPlanning.solve_quick(prob, x0, obs_polys)

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


direct_sol = PolyPlanning.solve_prob_direct_kkt(direct_prob, x0);

