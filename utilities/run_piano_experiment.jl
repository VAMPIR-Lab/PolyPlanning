using PolyPlanning
using JLD2
using Dates

# user options
is_saving = true
is_running_sep = true
is_running_kkt = false
is_loading_experiment = false
exp_file_date = "2024-05-30_1421"

# experiment parameters (ignored if is_loading_experiment)
n_maps = 2
n_x0s = 10
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
corridor_w_min = sqrt(2) * ego_length / 2
corridor_w_max = ego_length
pre_L_length_base = 3.0
init_x = pre_L_length_base + ego_length / 2
init_y_mean = -1.75
init_y_disturb_max = 1.0
init_θ_disturb_max = π / 2
data_dir = "data"
exp_name = "piano"
date_now = Dates.format(Dates.now(), "YYYY-mm-dd_HHMM")


if is_loading_experiment
    ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
else # generate ego_poly, x0s and maps
    @assert corridor_w_min^2 >= (ego_length / 2)^2 + ego_width^2
    @assert n_obs == 3

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
        init_x,
        init_y_mean,
        init_y_disturb_max
    )

    ego_poly = PolyPlanning.gen_ego_rect(; a=ego_width, b=ego_length)

    x0s = map(1:n_x0s) do i
        init_y = init_y_mean - init_y_disturb_max + 2 * init_y_disturb_max * rand()
        init_θ = -init_θ_disturb_max + 2 * init_θ_disturb_max * rand()
        [init_x, init_y, init_θ, 0, 0, 0]
    end

    maps = map(1:n_maps) do i
        width = corridor_w_min + (corridor_w_max - corridor_w_min) * rand() * 0
        # corridor entrance starts at the same x
        pre_L_length = pre_L_length_base - width / 2
        PolyPlanning.gen_L_corridor(; width, pre_L_length, post_L_length=1.0)
    end

    if is_saving
        jldsave("$data_dir/$(exp_name)_exp_$date_now.jld2"; ego_poly, x0s, maps, param)
    end
end

@info "Computing our solutions..."
start_t = time()
our_sols = PolyPlanning.multi_solve_ours(ego_poly, x0s, maps, param)
@info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
our_filt_by_success = PolyPlanning.filter_by_success(our_sols)
@info "our success rate $(length(our_filt_by_success.idx)/(n_maps*n_x0s)*100)%"

if is_saving
    jldsave("$data_dir/$(exp_name)_our_sols_$date_now.jld2"; our_sols)
end

if is_running_sep
    @info "Computing separating hyperplane solutions..."
    start_t = time()
    sep_sols = PolyPlanning.multi_solve_sep(ego_poly, x0s, maps, param)
    @info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
    sep_filt_by_success = PolyPlanning.filter_by_success(sep_sols)
    @info "sep success rate $(length(sep_filt_by_success.idx)/(n_maps*n_x0s)*100)%"

    if is_saving
        jldsave("$data_dir/$(exp_name)_sep_sols_$date_now.jld2"; sep_sols)
    end
end

if is_running_kkt
    @info "Computing direct KKT solutions..."
    start_t = time()
    kkt_sols = PolyPlanning.multi_solve_kkt(ego_poly, x0s, maps, param)
    @info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
    kkt_filt_by_success = PolyPlanning.filter_by_success(kkt_sols)
    @info "kkt success rate $(length(kkt_filt_by_success.idx)/(n_maps*n_x0s)*100)%"

    if is_saving
        jldsave("$data_dir/$(exp_name)_kkt_sols_$date_now.jld2"; kkt_sols)
    end
end

# load results from file
#our_sols, sep_sols, kkt_sols = PolyPlanning.load_results(exp_name, date_now; data_dir)

# visualize
PolyPlanning.visualize_multi(x0s, maps, our_sols, T, ego_poly; n_rows=2, n_cols=2, title_prefix = "ours")
PolyPlanning.visualize_multi(x0s, maps, sep_sols, T, ego_poly; n_rows=2, n_cols=2)

# easier table