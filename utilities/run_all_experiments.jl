using Dates
using PolyPlanning

# these turned out to be common settings in between experiments
is_saving = true
is_running_ours = true
is_running_sep = true
is_running_dcol = true
is_loading_exp = true # skip experiment generation and load from file
is_loading_res = true  # skip compute and load from file
exp_file_date = "2024-07-04_0428"
res_file_date = "2024-07-04_0428"
data_dir = "data"
date_now = Dates.format(Dates.now(), "YYYY-mm-dd_HHMM")
n_xu = 9
T = 20
dt = 0.2
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;
Qf = 2e-3 * PolyPlanning.I(2)
u1_max = 10.0
u2_max = 10.0
u3_max = π

@info date_now

@info "Running Simple packing experiment..."
include("run_simple_packing_exp.jl")

@info "Running simple gap experiment..."
include("run_simple_gap_exp.jl")

@info "Running Piano experiment..."
include("run_piano_exp.jl")

@info "Running Random packing experiment..."
include("run_random_packing_exp.jl")

@info "Running L Piano experiment..."
include("run_L_piano_exp.jl")

@info "Running Random L packing experiment..."
include("run_random_L_packing_exp.jl")