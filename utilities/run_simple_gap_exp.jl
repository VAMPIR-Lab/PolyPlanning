using PolyPlanning
using JLD2
using Dates
using GLMakie

# user options
is_saving = true
is_running_sep = true
is_running_dcol = true
is_running_kkt = false
is_loading_exp = false # skip experiment generation and load from file
is_loading_res = false  # skip compute and load from file
exp_file_date = "2024-05-30_2351"
res_file_date = "2024-05-30_2351"
exp_name = "simple_gap"
data_dir = "data"
date_now = Dates.format(Dates.now(), "YYYY-mm-dd_HHMM")

# experiment parameters (ignored if is_loading_exp or is_loading_res)
n_maps = 3 # number of maps
n_x0s = 100 # number of initial conditions
n_sides = 4 # 
n_obs = 2
n_xu = 9 # 6-state variable + control variable
T = 20 # timestep
dt = 0.2 #
Rf = 1e-3 * PolyPlanning.I(3); # penalty for control variable
Rf[3, 3] = Rf[3, 3] / 100.0;
Qf = 1e-2 * PolyPlanning.I(2) # penalty for translation
u1_max = 10.0
u2_max = 10.0
u3_max = π
n_sd_slots=4
init_x_mean = 5.0
init_y_mean = 0.0
init_x_disturb_max = 1.0
init_y_disturb_max = 1.0
ego_width = 0.5
ego_length = 2.0
gap_min = ego_width + 0.2
gap_max = ego_width * 2
gap_array = [gap_min, (gap_min + gap_max) / 2, gap_max]
gap_offset = 2.5

if is_loading_exp || is_loading_res
    ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
else # generate ego_poly, x0s and maps
    @assert init_x_mean - init_x_disturb_max - ego_length / 2 >= gap_offset
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
        n_sd_slots,
        init_x_mean,
        init_y_mean,
        init_x_disturb_max,
        init_y_disturb_max,
        data_dir,
        date_now,
        exp_name,
        ego_width,
        ego_length,
        gap_min,
        gap_max,
        gap_offset
    )

    # generate x0s and maps
    ego_poly = PolyPlanning.gen_ego_rect(; a=ego_width, b=ego_length)

    x0s = map(1:n_x0s) do i
        init_x = init_x_mean - init_x_disturb_max + 2 * init_x_disturb_max * rand()
        init_y = init_y_mean - init_y_disturb_max + 2 * init_y_disturb_max * rand()
        [init_x, init_y, -π + 2 * π * rand(), 0, 0, 0]
    end

    maps = map(1:n_maps) do i
        PolyPlanning.gen_gap(; width=gap_array[i], xs=gap_offset)
    end

    if is_saving
        exp_file_date = date_now
        jldsave("$data_dir/$(exp_name)_exp_$date_now.jld2"; ego_poly, x0s, maps, param)
    end
end

if is_loading_res
    our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.load_all(exp_name, exp_file_date, res_file_date; is_loading_sep=is_running_sep, is_loading_dcol=is_running_dcol, is_loading_kkt=is_running_kkt, data_dir)
else
    our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.compute_all(ego_poly, x0s, maps, param; is_saving, exp_name, date_now, exp_file_date, is_running_sep, is_running_dcol, is_running_kkt, data_dir)
end

# process
our_bins = PolyPlanning.process_into_bins(our_sols)
@info "$(length(our_bins.success.idx)/(param.n_maps*param.n_x0s)*100)% our success rate"

sep_bins = []
if is_running_sep
    sep_bins = PolyPlanning.process_into_bins(sep_sols)
    @info "$(length(sep_bins.success.idx)/(param.n_maps*param.n_x0s)*100)% sep success rate"
end

dcol_bins = []
if is_running_dcol
    dcol_bins = PolyPlanning.process_into_bins(dcol_sols)
    @info "$(length(dcol_bins.success.idx)/(param.n_maps*param.n_x0s)*100)% dcol success rate"
end

kkt_bins = []
if is_running_kkt
    kkt_bins = PolyPlanning.process_into_bins(kkt_sols)
    @info "$(length(kkt_bins.success.idx)/(param.n_maps*param.n_x0s)*100)% kkt success rate"
end

# tables
if is_running_sep
    PolyPlanning.print_stats(our_bins, sep_bins, param.n_maps, param.n_x0s; name="ours", ref_name="sep")
end

if is_running_dcol
    PolyPlanning.print_stats(our_bins, dcol_bins, param.n_maps, param.n_x0s; name="ours", ref_name="dcol")
end

if is_running_kkt
    PolyPlanning.print_stats(our_bins, kkt_bins, param.n_maps, param.n_x0s; name="ours", ref_name="kkt")
end
