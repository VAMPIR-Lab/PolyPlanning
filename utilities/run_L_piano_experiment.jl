using PolyPlanning
using JLD2
using Dates

# user options
is_saving = false
is_running_sep = false
is_running_kkt = false
is_loading_experiment = false
exp_file_date = "2024-05-30_1421"

# experiment parameters (ignored if is_loading_experiment)
n_maps = 1
n_x0s = 2
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
init_x = 5.0
init_y_max = 3.0
gap_min = 1.25
gap_max = 2.5
wall_xs = 3.0
data_dir = "data"
exp_name = "L_piano"
date_now = Dates.format(Dates.now(), "YYYY-mm-dd_HHMM")

if is_loading_experiment
    ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
else # generate ego_poly, x0s and maps
    @assert init_x >= wall_xs
    @assert n_obs == 2

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
        gap_min,
        gap_max,
        wall_xs,
        exp_name
    )

    # generate x0s and maps
    ego_poly = PolyPlanning.gen_ego_L()

    x0s = map(1:n_x0s) do i
        [init_x, -init_y_max + 2 * init_y_max * rand(), -π + 2 * π * rand(), 0, 0, 0]
    end

    maps = map(1:n_maps) do i
        PolyPlanning.gap_polys = PolyPlanning.gen_gap(; width=gap_min + (gap_max - gap_min) * rand(), xs=wall_xs)
    end

    if is_saving
        jldsave("$data_dir/$(exp_name)_exp_$date_now.jld2"; ego_poly, x0s, maps, param)
    end
end
# compute
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
#maps_idx = 1
#x0_idx = 1
#(fig, update_fig) = PolyPlanning.visualize_quick(x0s[x0_idx], T, ego_poly, maps[maps_idx])
#update_fig(our_sols[(maps_idx, x0_idx)].res.θ)
#display(fig)

#(fig, update_fig) = PolyPlanning.visualize_sep_planes(x0s[x0_idx], T, ego_poly, maps[maps_idx])
#update_fig(sep_sols[(maps_idx, x0_idx)].res.θ)
#display(fig)