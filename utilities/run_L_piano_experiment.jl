using PolyPlanning
using JLD2
using Dates

n_maps = 1
n_x0s = 2
n_sides = 4
n_obs = 2
n_xu = 9
T = 20
dt = 0.2
Rf = 1e-3 * PolyPlanning.I(3)
Rf[3, 3] = Rf[3, 3] / 10.0
Qf = 1e-2 * PolyPlanning.I(2)
u1_max = 1.0
u2_max = 1.0
u3_max = π
init_x = 5.0
init_y_max = 3.0
gap_min = 0.25
gap_max = 0.5
data_dir = "data"
exp_name = "packing"
date_now = PolyPlanning.Dates.format(PolyPlanning.Dates.now(), "YYYY-mm-dd_HHMM")

@assert init_x >= wall_w

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
    gap_min,
    gap_max
)

# generate x0s and maps
ego_poly = PolyPlanning.gen_ego_L();

x0s = map(1:n_x0s) do i
    [init_x, -init_y_max + 2 * init_y_max * rand(), -π + 2 * π * rand(), 0, 0, 0]
end

maps = map(1:n_maps) do i
    PolyPlanning.gap_polys = PolyPlanning.gen_gap(; length=gap_min + (gap_max - gap_min) * rand(), xs=3.0)
end

PolyPlanning.jldsave("$data_dir/$(exp_name)_exp_$date_now.jld2"; ego_poly, x0s, maps, param)

# compute
@info "Computing our solutions..."
start_t = time()
our_sols = PolyPlanning.multi_solve_ours(ego_poly, x0s, maps, param)
@info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
PolyPlanning.jldsave("$data_dir/$(exp_name)_our_sols_$date_now.jld2"; our_sols)

@info "Computing separating hyperplane solutions..."
start_t = time()
sep_sols = PolyPlanning.multi_solve_sep(ego_poly, x0s, maps, param)
@info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
PolyPlanning.jldsave("$data_dir/$(exp_name)_sep_sols_$date_now.jld2"; sep_sols)

@info "Computing direct KKT solutions..."
start_t = time()
kkt_sols = PolyPlanning.multi_solve_kkt(ego_poly, x0s, maps, param)
@info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
PolyPlanning.jldsave("$data_dir/$(exp_name)_kkt_sols_$date_now.jld2"; kkt_sols)

# load results from file
#ego_poly, x0s, maps, param, our_sols, sep_sols, kkt_sols = PolyPlanning.load_results(exp_name, date_now; data_dir)

# visualize
x0_idx = 1
maps_idx = 1
(fig, update_fig) = PolyPlanning.visualize_quick(x0s[x0_idx], T, ego_poly, maps[maps_idx])
update_fig(our_sols[(maps_idx, x0_idx)].res.θ)
display(fig)

# process results
#our_filt_by_success = PolyPlanning.filter_by_success(our_sols)
#our_v_sep = PolyPlanning.compute_sols_Δ(param.n_maps, param.n_x0s, our_sols, sep_sols)
