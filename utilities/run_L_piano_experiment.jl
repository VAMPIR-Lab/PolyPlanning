using PolyPlanning
using JLD2
using Dates

# user options
is_saving = true
is_running_sep = true
is_running_kkt = true
is_loading_exp = false # skip experiment generation and load from file
is_loading_res = false # skip compute and load from file
exp_file_date = "2024-05-30_2351"
res_file_date = "2024-05-30_2351"
exp_name = "L_piano"
data_dir = "data"
date_now = Dates.format(Dates.now(), "YYYY-mm-dd_HHMM")

# experiment parameters (ignored if is_loading_exp or is_loading_res)
n_maps = 3
n_x0s = 100
n_sides = 4
n_obs = 2
n_xu = 9
T = 20
dt = 0.2
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;
Qf = 5e-3 * PolyPlanning.I(2)
u1_max = 10.0
u2_max = 10.0
u3_max = π
init_x_mean = 5.0
init_y_mean = 0.0
init_x_disturb_max = 1.0
init_y_disturb_max = 3.0
gap_min = 1.25
gap_max = 2.5
gap_array = [gap_min, (gap_min + gap_max) / 2, gap_max]
wall_xs = 3.0

if is_loading_exp || is_loading_res
    ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
else # generate ego_poly, x0s and maps
    @assert init_x >= wall_xs
    @assert n_obs == 2
    @assert n_maps == length(gap_array)

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
        init_x_mean,
        init_y_mean,
        init_x_disturb_max,
        init_y_disturb_max,
        data_dir,
        date_now,
        gap_min,
        gap_max,
        wall_xs,
        exp_name
    )

    # generate x0s and maps
    ego_poly = PolyPlanning.gen_ego_L()

    x0s = map(1:n_x0s) do i
        init_x = init_x_mean - init_x_disturb_max + 2 * init_x_disturb_max * rand()
        init_y = init_y_mean - init_y_disturb_max + 2 * init_y_disturb_max * rand()
        [init_x, init_y, -π + 2 * π * rand(), 0, 0, 0]
    end

    maps = map(1:n_maps) do i
        PolyPlanning.gap_polys = PolyPlanning.gen_gap(; width=gap_min + (gap_max - gap_min) * rand(), xs=wall_xs)
    end

    if is_saving
        exp_file_date = date_now
        jldsave("$data_dir/$(exp_name)_exp_$date_now.jld2"; ego_poly, x0s, maps, param)
    end
end

if is_loading_res
    our_sols, sep_sols, kkt_sols = PolyPlanning.load_all(exp_name, res_file_date, exp_file_date; is_loading_sep=is_running_sep, is_loading_kkt=is_running_kkt, data_dir)
else
    our_sols, sep_sols, kkt_sols = PolyPlanning.compute_all(ego_poly, x0s, maps, param; is_saving, exp_name, date_now, exp_file_date, is_running_sep, is_running_kkt, data_dir)
end

# visualize
#PolyPlanning.visualize_multi(x0s, maps, our_sols, T, ego_poly; n_rows=3, n_cols=2, title_prefix="ours")
#PolyPlanning.visualize_multi(x0s, maps, sep_sols, T, ego_poly; n_rows=3, n_cols=2, title_prefix = "sep")
#PolyPlanning.visualize_multi(x0s, maps, kkt_sols, T, ego_poly; n_rows=3, n_cols=2, title_prefix = "kkt")

# process
