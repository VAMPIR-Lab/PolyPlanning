using PolyPlanning, JLD2, Statistics
data_dir = "data"
#experiment_file = jldopen("$(data_dir)/packing_experiment_2024-05-29_0603.jld2", "r")
#our_file = jldopen("$(data_dir)/packing_our_sols_2024-05-29_0603.jld2", "r")
#sep_file = jldopen("$(data_dir)/packing_sep_plane_sols_2024-05-29_0603.jld2", "r")
#kkt_file = jldopen("$(data_dir)/packing_direct_kkt_sols_2024-05-29_0603.jld2", "r")
experiment_file = jldopen("$(data_dir)/packing_experiment_2024-05-28_2140.jld2", "r")
our_file = jldopen("$(data_dir)/packing_our_sols_2024-05-28_2140.jld2", "r")
sep_file = jldopen("$(data_dir)/packing_sep_plane_sols_2024-05-28_2140.jld2", "r")
kkt_file = jldopen("$(data_dir)/packing_direct_kkt_sols_2024-05-28_2140.jld2", "r")
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

# simple processing
#our_sols_file = PolyPlanning.jldopen("$(data_dir)/packing_our_sols_2024-05-28_0007.jld2", "r")
#sep_sols_file = PolyPlanning.jldopen("$(data_dir)/packing_sep_plane_sols_2024-05-28_0007.jld2", "r")
#our_sols = our_sols_file["our_sols"]
#sep_plane_sols = sep_sols_file["sep_plane_sols"]

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

    if our_sol.mcp_success && sep_plane_sol.mcp_success# && direct_kkt_sol.mcp_success
        push!(idxs, idx)
        push!(our_times, our_sol.time)
        push!(our_x_dists, our_sol.x_dist)
        push!(sep_plane_times, sep_plane_sol.time)
        push!(sep_plane_x_dists, sep_plane_sol.x_dist)
        push!(direct_kkt_times, direct_kkt_sol.time)
        push!(direct_kkt_x_dists, direct_kkt_sol.x_dist)
    end
end

using Statistics

n_samples = length(idxs)

mean_our_times = mean(our_times);
mean_sep_times = mean(sep_plane_times);
#mean_kkt_times = mean(direct_kkt_times);

mean_our_times_CI = 1.96 * std(our_times) / sqrt(n_samples);
mean_sep_times_CI = 1.96 * std(sep_plane_times) / sqrt(n_samples);
#mean_kkt_times_CI = 1.96 * std(direct_kkt_times) / sqrt(n_samples);

@info "our times wrt direct kkt: $(round(mean_our_times/mean_sep_times*100; sigdigits=3))% (±$(round(mean_our_times_CI/mean_sep_times*100; sigdigits=3)))"

@info "sep plane times wrt direct kkt: $(round(mean_sep_times/mean_sep_times*100; sigdigits=3))% (±$(round(mean_sep_times_CI/mean_sep_times*100; sigdigits=3)))"

mean_our_x_dists = mean(our_x_dists);
mean_sep_x_dists = mean(sep_plane_x_dists);
#mean_kkt_x_dists = mean(direct_kkt_x_dists);

mean_our_x_dists_CI = 1.96 * std(our_x_dists) / sqrt(n_samples);
mean_sep_x_dists_CI = 1.96 * std(sep_plane_x_dists) / sqrt(n_samples);
#mean_kkt_x_dists_CI = 1.96 * std(direct_kkt_x_dists) / sqrt(n_samples);

@info "our x dists wrt direct kkt: $(round(mean_our_x_dists/mean_sep_x_dists*100; sigdigits=3))% (±$(round(mean_our_x_dists_CI/mean_sep_x_dists*100; sigdigits=3)))"

@info "sep plane x dists wrt direct kkt: $(round(mean_sep_x_dists/mean_sep_x_dists*100; sigdigits=3))% (±$(round(mean_sep_x_dists_CI/mean_sep_x_dists*100; sigdigits=3)))"
