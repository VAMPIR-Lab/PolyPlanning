using PolyPlanning
using JLD2
using Dates
using Statistics

# user options
is_saving = true
is_running_sep = true
is_running_kkt = false
is_loading_exp = false # skip experiment generation and load from file
is_loading_res = false # skip compute and load from file
exp_file_date = "2024-05-31_0108"
res_file_date = "2024-05-31_0108"
exp_name = "piano"
data_dir = "data"
date_now = Dates.format(Dates.now(), "YYYY-mm-dd_HHMM")

# experiment parameters (ignored if is_loading_exp or is_loading_res)
n_maps = 3
n_x0s = 8
n_sides = 4
n_obs = 3
n_xu = 9
T = 20
dt = 0.2
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;
Qf = 5e-3 * PolyPlanning.I(2)
u1_max = 10.0
u2_max = 10.0
u3_max = π
ego_width = 0.5
ego_length = 2.0
corridor_w_min = sqrt((ego_length / 2)^2 + ego_width^2)
corridor_w_max = ego_length
corridor_w_array = [corridor_w_min, (corridor_w_min + corridor_w_max) / 2, corridor_w_max]
pre_L_length_base = 4.0
post_L_length = 1.0
init_x_min = pre_L_length_base + ego_length / 2
init_y_mean = -post_L_length - corridor_w_min / 2
init_x_disturb_max = corridor_w_min / 2
init_y_disturb_max = corridor_w_min / 2
init_θ_disturb_max = π / 4

if is_loading_exp || is_loading_res
    ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
else # generate ego_poly, x0s and maps
    @assert corridor_w_min^2 >= (ego_length / 2)^2 + ego_width^2
    @assert n_obs == 3
    @assert n_maps == length(corridor_w_array)

    param = (;
        n_maps,
        n_x0s,
        n_sides,
        n_obs,
        n_xu,
        T,
        dt,
        Rf,
        Qf,
        u1_max,
        u2_max,
        u3_max,
        data_dir,
        date_now,
        exp_name,
        ego_width,
        ego_length,
        corridor_w_min,
        corridor_w_max,
        init_x_min,
        init_y_mean,
        init_x_disturb_max,
        init_y_disturb_max,
        init_θ_disturb_max
    )

    ego_poly = PolyPlanning.gen_ego_rect(; a=ego_width, b=ego_length)

    x0s = map(1:n_x0s) do i
        init_x = init_x_min + init_x_disturb_max * rand()
        init_y = init_y_mean - init_y_disturb_max + 2 * init_y_disturb_max * rand()
        init_θ = -init_θ_disturb_max + 2 * init_θ_disturb_max * rand()
        [init_x, init_y, init_θ, 0, 0, 0]
    end

    maps = map(1:n_maps) do i
        #width = corridor_w_min + (corridor_w_max - corridor_w_min) * rand()
        width = corridor_w_array[i]
        # corridor entrance starts at the same x
        pre_L_length = pre_L_length_base - width / 2
        PolyPlanning.gen_L_corridor(; width, pre_L_length, post_L_length)
    end

    if is_saving
        exp_file_date = date_now
        jldsave("$data_dir/$(exp_name)_exp_$exp_file_date.jld2"; ego_poly, x0s, maps, param)
    end
end

if is_loading_res
    our_sols, sep_sols, kkt_sols = PolyPlanning.load_all(exp_name, res_file_date, exp_file_date; is_loading_sep=is_running_sep, is_loading_kkt=is_running_kkt, data_dir)
else
    our_sols, sep_sols, kkt_sols = PolyPlanning.compute_all(ego_poly, x0s, maps, param; is_saving, exp_name, date_now, exp_file_date, is_running_sep, is_running_kkt, data_dir)
end

# visualize
PolyPlanning.visualize_multi_quick(x0s, maps, our_sols, T, ego_poly; n_rows=3, n_cols=8, title_prefix="ours")
#PolyPlanning.visualize_multi(x0s, maps, sep_sols, T, ego_poly; n_rows=3, n_cols=2, title_prefix = "sep")
#PolyPlanning.visualize_multi(x0s, maps, kkt_sols, T, ego_poly; n_rows=3, n_cols=2, title_prefix = "kkt")

# process
our_only_mcp_success = PolyPlanning.filter_by_mcp_success(our_sols)
our_only_task_success = PolyPlanning.filter_by_task_success(our_sols; task_radius=0.5)
sep_only_mcp_success = PolyPlanning.filter_by_mcp_success(sep_sols)
sep_only_task_success = PolyPlanning.filter_by_task_success(sep_sols; task_radius=0.5)

our_mcp_success_ratio = length(our_only_mcp_success.idx) / (param.n_maps * param.n_x0s)
our_task_success_ratio = length(our_only_task_success.idx) / (param.n_maps * param.n_x0s)
sep_mcp_success_ratio = length(sep_only_mcp_success.idx) / (param.n_maps * param.n_x0s)
sep_task_success_ratio = length(sep_only_task_success.idx) / (param.n_maps * param.n_x0s)

our_Δ_wrt_sep_mcp = PolyPlanning.compute_sols_Δ_mcp(param.n_maps, param.n_x0s, our_sols, sep_sols)
our_Δ_wrt_sep_task = PolyPlanning.compute_sols_Δ_task(param.n_maps, param.n_x0s, our_sols, sep_sols)

mean_x_dist_Δ_mcp = mean(our_Δ_wrt_sep_mcp.x_dist_Δ)
mean_time_Δ_mcp = mean(our_Δ_wrt_sep_mcp.time_Δ)
mean_x_dist_Δ_task = mean(our_Δ_wrt_sep_task.x_dist_Δ)
mean_time_Δ_task = mean(our_Δ_wrt_sep_task.time_Δ)

println("                    ours     sep ")
println("only mcp success:")
println("mcp success  $(round(our_mcp_success_ratio*100; sigdigits=2))%  $(round(sep_mcp_success_ratio*100; sigdigits=2))%")
println("avg Δ x dist  $(round(mean_x_dist_Δ_mcp; sigdigits=2))  0.0")
println("avg Δ time  $(round(mean_time_Δ_mcp; sigdigits=2))  0.0")
println("only task success:")
println("task success  $(round(our_task_success_ratio*100; sigdigits=2))%  $(round(sep_task_success_ratio*100; sigdigits=2))%")
println("avg Δ x dist $(round(mean_x_dist_Δ_task; sigdigits=2))  0.0")
println("avg Δ time  $(round(1mean_time_Δ_task; sigdigits=2))  0.0")