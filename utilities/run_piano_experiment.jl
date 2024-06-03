using PolyPlanning
using JLD2
using Dates
using Statistics

# user options
is_saving = true
is_running_sep = true
is_running_kkt = true
is_loading_exp = false # skip experiment generation and load from file
is_loading_res = false # skip compute and load from file
exp_file_date = "2024-06-03_0140"
res_file_date = "2024-06-03_0140"
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
pre_L_length_base = 5.0
post_L_length = 3.0
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

# process
our_bins = PolyPlanning.process_into_bins(our_sols)
@info "$(length(our_bins.successes.idx)/(param.n_maps*param.n_x0s)*100)% nonsmooth success rate"

sep_bins = []
if is_running_sep
    sep_bins = PolyPlanning.process_into_bins(sep_sols)
    @info "$(length(sep_bins.successes.idx)/(param.n_maps*param.n_x0s)*100)% sep success rate"
end

kkt_bins = []
if is_running_kkt
    kkt_bins = PolyPlanning.process_into_bins(kkt_sols)
    @info "$(length(kkt_bins.successes.idx)/(param.n_maps*param.n_x0s)*100)% kkt success rate"
end

# visualize
#PolyPlanning.visualize_multi(x0s, maps, our_sols, T, ego_poly; n_rows=3, n_cols=8, type="nonsmooth")
#PolyPlanning.visualize_multi(x0s, maps, sep_sols, T, ego_poly; n_rows=3, n_cols=8, type="sep_planes")
#PolyPlanning.visualize_multi(x0s, maps, kkt_sols, T, ego_poly; n_rows=3, n_cols=2, type = "direct_kkt")

#PolyPlanning.visualize_multi(x0s, maps, our_sols, our_bins.successes, T, ego_poly; n_rows=3, n_cols=8, type="nonsmooth")
#PolyPlanning.visualize_multi(x0s, maps, sep_sols, sep_bins.successes, T, ego_poly; n_rows=3, n_cols=8, type="sep_planes")
#PolyPlanning.visualize_multi(x0s, maps, kkt_sols, kkt_bins.successes, T, ego_poly; n_rows=3, n_cols=8, type="direct_kkt")

# tables
PolyPlanning.print_table(our_bins, param.n_maps, param.n_x0s; title="our")

if is_running_sep
    PolyPlanning.print_table(sep_bins, param.n_maps, param.n_x0s, title="sep")
end

if is_running_kkt
    PolyPlanning.print_table(kkt_bins, param.n_maps, param.n_x0s, title="kkt")
end

if is_running_sep
    our_v_sep_locals = PolyPlanning.compute_Δ_time_cost(our_bins.successes, sep_bins.successes)
    our_v_sep_fails = PolyPlanning.compute_Δ_time_cost(our_bins.fails, sep_bins.fails)
    our_v_sep_bin = (successes=our_v_sep_locals, fails=our_v_sep_fails)

    PolyPlanning.print_table(our_v_sep_bin, param.n_maps, param.n_x0s, title="our v sep")

    if is_running_kkt
        our_v_kkt_locals = PolyPlanning.compute_Δ_time_cost(kkt_bins.successes, sep_bins.successes)
        our_v_kkt_fails = PolyPlanning.compute_Δ_time_cost(kkt_bins.fails, sep_bins.fails)
        our_v_kkt_bin = (successes=our_v_kkt_locals, fails=our_v_kkt_fails)

        PolyPlanning.print_table(our_v_kkt_bin, param.n_maps, param.n_x0s, title="our v kkt")
    end
end
