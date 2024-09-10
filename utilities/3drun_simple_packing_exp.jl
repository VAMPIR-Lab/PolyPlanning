using PolyPlanning
using JLD2
using Dates
using GLMakie
using LinearAlgebra

# user options
is_saving = true
is_running_ourss = true
is_running_sep = true
is_running_dcol = false
is_running_kkt = false
is_loading_exp = true # skip experiment generation and load from file
is_loading_res = true  # skip compute and load from file
exp_file_date = "2024-08-19_0033"
res_file_date = "2024-08-19_0033"
exp_name = "simple_packing"
data_dir = "data"
date_now = Dates.format(Dates.now(), "YYYY-mm-dd_HHMM")

# experiment parameters (ignored if is_loading_exp or is_loading_res)
n_maps = 10 # number of maps
n_x0s = 100 # number of initial conditions
n_sides = 4 # 
n_obs = 1
n_xu = 12 # 6-state variable + control variable
T = 2 # timestep
dt = 2. #
R_cost = 1e-3 * PolyPlanning.I(6) # penalty for control variable
R_cost[4:6, 4:6] = R_cost[4:6, 4:6] / 100.0
Q_cost=2e-3 * PolyPlanning.I(3) # penalty for translation
u1_max=10.0
u2_max=10.0
u3_max=10.0
u4_max=π
u5_max=π
u6_max=π

r=rand
if is_loading_exp || is_loading_res
    ego_poly, x0s, maps, params = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
    param = (;
        n_maps,
        n_x0s,
        n_sides,
        n_obs,
        n_xu,
        T,
        dt,
        R_cost,
        Q_cost,
        u1_max=10.0,
        u2_max=10.0,
        u3_max=10.0,
        u4_max=π,
        u5_max=π,
        u6_max=π,
        data_dir,
        date_now,
        exp_name
    )
else # generate ego_poly, x0s and maps

    param = (;
        n_maps,
        n_x0s,
        n_sides,
        n_obs,
        n_xu,
        T,
        dt,
        R_cost,
        Q_cost,
        u1_max=10.0,
        u2_max=10.0,
        u3_max=10.0,
        u4_max=π,
        u5_max=π,
        u6_max=π,
        data_dir,
        date_now,
        exp_name
    )

    # generate x0s and maps
    Ve = [[-1-r(), -1-r(), -1-r()], [1+r(), -1-r(), -1-r()], [0, 1+r(), -1-r()], [0, 0, 1+r()]]
    Pe = PolyPlanning.ConvexPolygon3D(Ve)
    ego_poly = [Pe]

    x0s = map(1:n_x0s) do i
        axis = rand(3)
        mrp = axis / norm(axis) * (0.8+0.4*r()) # norm(mrp)∈[0.8, 1.2]
        # mrp = (ones(3)*0.8+r(3)*0.4)/sqrt(3) # tan(θ/4)∈[0.8, 1.2] θ∈[2.70, 3.50]
        trans =zeros(3) + [5,2,3] + r(3)
        [trans; mrp; zeros(6)]
    end


    maps = map(1:n_maps) do i
        Vo = [[-1-r(), -1-r(), -1-r()], [1+r(), -1-r(), -1-r()], [0, 1+r(), -1-r()], [0, 0, 5r()]]
        Po = PolyPlanning.ConvexPolygon3D(Vo)
        [Po]
    end

    if is_saving
        exp_file_date = date_now
        jldsave("$data_dir/$(exp_name)_exp_$date_now.jld2"; ego_poly, x0s, maps, param)
    end
end

if is_loading_res
    our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.load_all(exp_name, exp_file_date, res_file_date; is_loading_our=is_running_ourss, is_loading_sep=is_running_sep, is_loading_dcol=is_running_dcol, is_loading_kkt=is_running_kkt, data_dir)
else
    our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.compute_all_3d(ego_poly, x0s, maps, param; is_saving, exp_name, date_now, exp_file_date,is_running_ourss, is_running_sep, is_running_dcol, is_running_kkt, data_dir)
end


# process
our_bins = []
if is_running_ours
    our_bins = PolyPlanning.process_into_bins(our_sols)
    #@info "$(length(our_bins.success.idx)/(length(our_sols))*100)% our success rate"
end

sep_bins = []
if is_running_sep
    sep_bins = PolyPlanning.process_into_bins(sep_sols)
    #@info "$(length(sep_bins.success.idx)/(length(sep_sols))*100)% sep success rate"
end

dcol_bins = []
if is_running_dcol
    dcol_bins = PolyPlanning.process_into_bins(dcol_sols)
    #@info "$(length(dcol_bins.success.idx)/(length(dcol_sols))*100)% dcol success rate"
end

kkt_bins = []
if is_running_kkt
    kkt_bins = PolyPlanning.process_into_bins(kkt_sols)
    #@info "$(length(kkt_bins.success.idx)/(length(kkt_sols))*100)% kkt success rate"
end



# tables
if is_running_ours
    PolyPlanning.print_stats(our_bins, length(our_sols); name="ours")
end

if is_running_sep
    PolyPlanning.print_stats(sep_bins, length(sep_sols); name="sep")
end

if is_running_dcol
    PolyPlanning.print_stats(dcol_bins, length(dcol_sols); name="dcol")
end

if is_running_kkt
    PolyPlanning.print_stats(kkt_bins, length(kkt_sols); name="kkt")
end

if is_running_ours
    if is_running_sep
        PolyPlanning.print_stats(our_bins, sep_bins, length(our_sols); name="ours", ref_name="sep")
    end

    if is_running_dcol
        PolyPlanning.print_stats(our_bins, dcol_bins, length(our_sols); name="ours", ref_name="dcol")
    end

    if is_running_kkt
        PolyPlanning.print_stats(our_bins, kkt_bins, length(our_sols); name="ours", ref_name="kkt")
    end
end