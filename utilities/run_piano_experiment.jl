using PolyPlanning
using JLD2
using Dates

n_maps = 1
n_x0s = 1
n_sides = 4
n_obs = 3
n_xu = 9
T = 20
dt = 0.2
Rf = 1e-3 * PolyPlanning.I(3)
Rf[3, 3] = Rf[3, 3] / 10.0
Qf = 1e-2 * PolyPlanning.I(2)
u1_max = 1.0
u2_max = 1.0
u3_max = π
init_x = 2.0
init_y_shift = -1.5
init_y_max = .5
#wall_w = 5.0
#wall_l = 5.0
data_dir = "data"
exp_name = "piano"
date_now = Dates.format(Dates.now(), "YYYY-mm-dd_HHMM")

#@assert init_x >= wall_w
@assert n_obs == 3

param = (; n_maps,
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
    wall_w,
    wall_l,
    exp_name
)

# generate x0s and maps
ego_poly = PolyPlanning.gen_ego_rect(; l_multip=2.0);

x0s = map(1:n_x0s) do i
    [init_x, init_y_shift - init_y_max + 2 * init_y_max * rand(), -π + 2 * π * rand(), 0, 0, 0]
end

maps = map(1:n_maps) do i
    PolyPlanning.gen_L_corridor(; width=1.0, pre_L_length=1.0, post_L_length=1.0)
end

#jldsave("$data_dir/$(exp_name)_exp_$date_now.jld2"; ego_poly, x0s, maps, param)

# compute
@info "Computing our solutions..."
start_t = time()
our_sols = PolyPlanning.multi_solve_ours(ego_poly, x0s, maps, param)
@info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
#jldsave("$data_dir/$(exp_name)_our_sols_$date_now.jld2"; our_sols)

#@info "Computing separating hyperplane solutions..."
#start_t = time()
#sep_sols = PolyPlanning.multi_solve_sep(ego_poly, x0s, maps, param)
#@info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
#jldsave("$data_dir/$(exp_name)_sep_sols_$date_now.jld2"; sep_sols)

#@info "Computing direct KKT solutions..."
#start_t = time()
#kkt_sols = PolyPlanning.multi_solve_kkt(ego_poly, x0s, maps, param)
#@info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
#jldsave("$data_dir/$(exp_name)_kkt_sols_$date_now.jld2"; kkt_sols)

## load results from file
#ego_poly, x0s, maps, param, our_sols, sep_sols, kkt_sols = PolyPlanning.load_results(exp_name, date_now; data_dir)

# process results
#our_filt_by_success = PolyPlanning.filter_by_success(our_sols)
#our_v_sep = PolyPlanning.compute_sols_Δ(param.n_maps, param.n_x0s, our_sols, sep_sols)

# visualize
#maps_idx = 1
#x0_idx = 3
#(fig, update_fig) = PolyPlanning.visualize_quick(x0s[x0_idx], T, ego_poly, maps[maps_idx])
#update_fig(our_sols[(maps_idx, x0_idx)].res.θ)
#display(fig)

#(fig, update_fig) = PolyPlanning.visualize_sep_planes(x0s[x0_idx], T, ego_poly, maps[maps_idx])
#update_fig(sep_sols[(maps_idx, x0_idx)].res.θ)
#display(fig)


