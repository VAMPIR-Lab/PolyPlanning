using Dates
using PolyPlanning

# these turned out to be common settings in between experiments
is_saving = true
is_running_sep = true
is_running_dcol = true
is_loading_exp = false # skip experiment generation and load from file
is_loading_res = false  # skip compute and load from file
exp_file_date = "2024-06-21_1348"
res_file_date = "2024-06-17_1454"
exp_name = "simple_packing"
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

include("run_simple_gap_exp.jl")
include("run_simple_packing_exp.jl")
include("run_piano_exp.jl")
include("run_random_packing_exp.jl")
include("run_L_through_gap_exp.jl")
include("run_random_L_packing_exp.jl")