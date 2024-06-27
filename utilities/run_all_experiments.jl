@info "simple packing experiment"
include("run_simple_packing_exp.jl")

@info "simple gap experiment"
include("run_simple_gap_exp.jl")

@info "piano experiment"
include("run_piano_exp.jl")

@info "random packing experiment"
include("run_random_packing_exp.jl")

@info "L through gap experiment"
include("run_L_through_gap_exp.jl")

@info "random L packing experiment"
include("run_random_L_packing_exp.jl")