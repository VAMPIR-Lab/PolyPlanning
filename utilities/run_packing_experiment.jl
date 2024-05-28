# random initial conditions
# random walls
# random ego shape? (later)

# metrics: 
# time to solve
# success from mcp status
# success (from constraints) (later)
# x dist

using PolyPlanning
using ProgressMeter

n_walls = 1
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

# gen experiment data
ego_rect = PolyPlanning.gen_ego_rect(; l_multip=2.0);

#walls = map(1:n_walls) do i
#    PolyPlanning.gen_packing_wall(n_obs, n_sides; w=1.5, l=2.0, max_overlap=0.5)
#end
#PolyPlanning.plot_polys(walls[1])

x0s = map(1:n_x0s) do i
    [init_x, -init_y_max + 2 * init_y_max * rand(), -π + 2 * π * rand(), 0, 0, 0]
end

PolyPlanning.jldsave("$(data_dir)/packing_experiment_$(date_now).jld2"; ego_rect, walls, x0s, params)

# ours
our_prob = PolyPlanning.setup_quick(
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
p = Progress(n_walls * n_x0s, dt=1.0)

for (i, wall) in enumerate(walls)
    for (j, x0) in enumerate(x0s)
        res = PolyPlanning.solve_quick(our_prob, x0, wall; is_displaying=false)
        mcp_success = res.status == PolyPlanning.PATHSolver.MCP_Solved
        time = res.info.total_time
        x_dist = res.z[(T-1)*n_xu+1]

        our_sols[i, j] = (; mcp_success, time, x_dist, res)
        next!(p)
    end
end
PolyPlanning.jldsave("$(data_dir)/packing_our_sols_$(date_now).jld2"; our_sols)

@info "computing sep plane sols"
sep_plane_sols = Dict()
p = Progress(n_walls * n_x0s, dt=1.0)

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
        next!(p)
    end
end
PolyPlanning.jldsave("$(data_dir)/packing_sep_plane_sols_$(date_now).jld2"; sep_plane_sols)


@info "computing direct kkt sols"
direct_kkt_sols = Dict()
p = Progress(n_walls * n_x0s, dt=1.0)

for (i, wall) in enumerate(walls)
    direct_kkt_prob = PolyPlanning.setup_direct_kkt(
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
        res = PolyPlanning.solve_prob_direct_kkt(direct_kkt_prob, x0; is_displaying=false)
        mcp_success = res.status == PolyPlanning.PATHSolver.MCP_Solved
        time = res.info.total_time
        x_dist = res.z[(T-1)*n_xu+1]

        direct_kkt_sols[i, j] = (; mcp_success, time, x_dist, res)
        next!(p)
    end
end

PolyPlanning.jldsave("$(data_dir)/packing_direct_kkt_sols_$(date_now).jld2"; direct_kkt_sols)

# simple processing
idxs = []
our_times = []
our_x_dists = []
our_idxs = []
sep_plane_times = []
sep_plane_x_dists = []
direct_kkt_times = []
direct_kkt_x_dists = []

for (idx, our_sol) in our_sols
    sep_plane_sol = sep_plane_sols[idx]
    direct_kkt_sol = direct_kkt_sols[idx]

    if our_sol.mcp_success && sep_plane_sol.mcp_success && direct_kkt_sol.mcp_success
        push!(idxs, idx)
        push!(our_times, our_sol.time)
        push!(our_x_dists, our_sol.x_dist)
        push!(sep_plane_times, sep_plane_sol.time)
        push!(sep_plane_x_dists, sep_plane_sol.x_dist)
        push!(direct_kkt_times, direct_kkt_sol.time)
        push!(direct_kkt_x_dists, direct_kkt_sol.x_dist)
    end
end

Main.@infiltrate
i_wall = 1
j_x0 = 2
(fig, update_fig) = PolyPlanning.visualize_quick(x0s[j_x0], T, ego_rect, walls[i_wall])
update_fig(our_sols[i_wall, j_x0].res.θ)
display(fig)

i_wall = 1
j_x0 = 2
(fig, update_fig) = PolyPlanning.visualize_sep_planes(x0s[j_x0], T, ego_rect, walls[i_wall])
update_fig(sep_plane_sols[i_wall, j_x0].res.θ)
display(fig)

i_wall = 1
j_x0 = 2
(fig, update_fig) = PolyPlanning.visualize_direct_kkt(x0s[j_x0], T, ego_rect, walls[i_wall])
update_fig(direct_kkt_sols[i_wall, j_x0].res.θ)
display(fig)
