using PolyPlanning
using JLD2
using Dates
using GLMakie

# user options
is_saving = false
is_running_sep = false
is_running_kkt = false
is_loading_exp = true # skip experiment generation and load from file
is_loading_res = false  # skip compute and load from file
exp_file_date = "2024-06-11_0929"
res_file_date = "2024-06-11_0929"
exp_name = "L simple" # L simple or simple
data_dir = "data"
date_now = Dates.format(Dates.now(), "YYYY-mm-dd_HHMM")

# experiment parameters (ignored if is_loading_exp or is_loading_res)
n_maps = 3 # number of maps
n_x0s = 12 # number of initial conditions
n_sides = 4 # 
n_obs = 1
n_xu = 9 # 6-state variable + control variable
T = 20 # timestep
dt = 0.2 #
Rf = 1e-3 * PolyPlanning.I(3); # penality for control variable
Rf[3, 3] = Rf[3, 3] / 100.0;
Qf = 1e-3 * PolyPlanning.I(2) # penality for translation
u1_max = 10.0
u2_max = 10.0
u3_max = π
init_x = 3.0
init_y_max = 1.0
ego_width = 0.5
ego_length = 1.0

if is_loading_exp || is_loading_res
    ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
else # generate ego_poly, x0s and maps
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
        init_x,
        init_y_max,
        data_dir,
        date_now,
        exp_name,
        ego_width,
        ego_length
    )
    # generate x0s and maps
    # ego_poly = PolyPlanning.gen_ego_rect(; a=ego_width, b=ego_length)
    ego_poly = PolyPlanning.PolyPlanning.gen_ego_L()

    x0s = map(1:n_x0s) do i
        [init_x / 2 * (1 + 2 * rand()), -init_y_max + 2 * init_y_max * rand(), -π + 2 * π * rand(), 0, 0, 0]
    end

    maps = map(1:n_maps) do i
        PolyPlanning.gen_rect_obs(; a=0.2 + 0.1*rand())
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
PolyPlanning.visualize_multi(x0s, maps, our_sols, T, ego_poly; n_rows=3, n_cols=4, title_prefix="ours")

#PolyPlanning.visualize_multi(x0s, maps, sep_sols, T, ego_poly; n_rows=3, n_cols=2, title_prefix = "sep")
#PolyPlanning.visualize_multi(x0s, maps, kkt_sols, T, ego_poly; n_rows=3, n_cols=2, title_prefix = "kkt")

# process