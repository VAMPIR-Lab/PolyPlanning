using PolyPlanning, JLD2, Statistics
data_dir = "data"
experiment_file = jldopen("$(data_dir)/packing_experiment_2024-05-29_0603.jld2", "r")
our_file = jldopen("$(data_dir)/packing_our_sols_2024-05-29_0603.jld2", "r")
sep_file = jldopen("$(data_dir)/packing_sep_plane_sols_2024-05-29_0603.jld2", "r")
kkt_file = jldopen("$(data_dir)/packing_direct_kkt_sols_2024-05-29_0603.jld2", "r")
ego_rect = experiment_file["ego_rect"]
walls = experiment_file["walls"]
x0s = experiment_file["x0s"]
params = experiment_file["params"]
our_sols = our_file["our_sols"]
sep_sols = sep_file["sep_plane_sols"]
kkt_sols = kkt_file["direct_kkt_sols"]


our_v_sep_indexes = []
our_v_sep_x_dist = []
our_v_sep_time = []

for (i, w) in enumerate(walls)
    for (j, x0) in enumerate(x0s)
        our_sol = our_sols[(i, j)]
        sep_sol = sep_sols[(i, j)]
		kkt_sol = kkt_sols[(i, j)]
        if our_sol.mcp_success && sep_sol.mcp_success
			push!(our_v_sep_indexes, (i, j))
            push!(our_v_sep_x_dist, our_sol.x_dist - sep_sol.x_dist)
			push!(our_v_sep_time, our_sol.time - sep_sol.time)
        end
    end
end
kkt_v_sep_indexes = []
kkt_v_sep_x_dist = []
kkt_v_sep_time = []
for (i, w) in enumerate(walls)
    for (j, x0) in enumerate(x0s)
        sep_sol = sep_sols[(i, j)]
		kkt_sol = kkt_sols[(i, j)]
        if sep_sol.mcp_success && kkt_sol.mcp_success
			push!(kkt_v_sep_indexes, (i, j))
			push!(kkt_v_sep_x_dist, kkt_sol.x_dist - sep_sol.x_dist)
			push!(kkt_v_sep_time, kkt_sol.time - sep_sol.time)
        end
    end
end


our_v_sep_samples = length(our_v_sep_indexes)
our_v_sep_x_dist_CI = 1.96 * std(our_v_sep_x_dist) / sqrt(our_v_sep_samples);
our_v_sep_time_CI = 1.96 * std(our_v_sep_time) / sqrt(our_v_sep_samples);
kkt_v_sep_x_dist_CI = 1.96 * std(kkt_v_sep_x_dist) / sqrt(our_v_sep_samples);
kkt_v_sep_time_CI = 1.96 * std(kkt_v_sep_time) / sqrt(our_v_sep_samples);