function multi_solve_ours(ego_poly, x0s, maps, param)
    @argcheck length(x0s) == param.n_x0s
    @argcheck length(maps) == param.n_maps
    for map in maps
        @argcheck length(map) == param.n_obs
    end

    prob = setup_quick(
        ego_poly;
        param.T,
        param.dt,
        Q=0.0 * [1.0 0; 0 1], # disabled final cost
        q=[0, 0.0],
        param.Rf,
        param.Qf,
        param.u1_max,
        param.u2_max,
        param.u3_max,
        sides_per_poly=4,
        derivs_per_sd=4,
        derivs_per_fv=4,
        param.n_obs
    )

    sols = Dict()
    p = Progress(param.n_maps * param.n_x0s, dt=1.0)
    for (i, map) in enumerate(maps)
        for (j, x0) in enumerate(x0s)
            res = solve_quick(prob, x0, map; is_displaying=false)
            mcp_success = res.status == PATHSolver.MCP_Solved
            time = res.info.total_time
            x_dist = res.z[(param.T-1)*param.n_xu+1]
            sols[i, j] = (; mcp_success, time, x_dist, res)
            next!(p)
        end
    end
    sols
end

function multi_solve_sep(ego_poly, x0s, maps, param)
    @argcheck length(x0s) == param.n_x0s
    @argcheck length(maps) == param.n_maps
    for map in maps
        @argcheck length(map) == param.n_obs
    end

    sols = Dict()
    p = Progress(param.n_maps * param.n_x0s, dt=1.0)
    for (i, map) in enumerate(maps)
        prob = setup_sep_planes(
            ego_poly,
            map;
            param.T,
            param.dt,
            param.Rf,
            param.Qf,
            param.u1_max,
            param.u2_max,
            param.u3_max
        )

        for (j, x0) in enumerate(x0s)
            res = solve_prob_sep_planes(prob, x0; is_displaying=false)
            mcp_success = res.status == PATHSolver.MCP_Solved
            time = res.info.total_time
            x_dist = res.z[(param.T-1)*param.n_xu+1]
            sols[i, j] = (; mcp_success, time, x_dist, res)
            next!(p)
        end
    end
    sols
end

function multi_solve_kkt(ego_poly, x0s, maps, param)
    @argcheck length(x0s) == param.n_x0s
    @argcheck length(maps) == param.n_maps
    for map in maps
        @argcheck length(map) == param.n_obs
    end

    sols = Dict()
    p = Progress(param.n_maps * param.n_x0s, dt=1.0)
    for (i, map) in enumerate(maps)
        prob = setup_direct_kkt(
            ego_poly,
            map;
            param.T,
            param.dt,
            param.Rf,
            param.Qf,
            param.u1_max,
            param.u2_max,
            param.u3_max
        )

        for (j, x0) in enumerate(x0s)
            res = solve_prob_direct_kkt(prob, x0; is_displaying=false)
            mcp_success = res.status == PATHSolver.MCP_Solved
            time = res.info.total_time
            x_dist = res.z[(param.T-1)*param.n_xu+1]
            sols[i, j] = (; mcp_success, time, x_dist, res)
            next!(p)
        end
    end
    sols
end

function load_results(name, date; data_dir="data")
    exp_file = jldopen("$data_dir/$(name)_exp_$date.jld2", "r")
    our_file = jldopen("$data_dir/$(name)_our_sols_$date.jld2", "r")
    sep_file = jldopen("$data_dir/$(name)_sep_sols_$date.jld2", "r")
    kkt_file = jldopen("$data_dir/$(name)_kkt_sols_$date.jld2", "r")

    ego_poly = exp_file["ego_poly"]
    x0s = exp_file["x0s"]
    maps = exp_file["maps"]
    param = exp_file["param"]
    our_sols = our_file["our_sols"]
    sep_sols = sep_file["sep_sols"]
    kkt_sols = kkt_file["kkt_sols"]

    (; ego_poly, x0s, maps, param, our_sols, sep_sols, kkt_sols)
end

function filter_by_success(sols)
    idx = []
    times = []
    x_dists = []

    for (i, sol) in sols
        if sol.mcp_success
            push!(idx, i)
            push!(times, sol.time)
            push!(x_dists, sol.x_dist)
        end
    end
    (; idx, times, x_dists)
end

function compute_sols_Δ(n_maps, n_x0s, sols, ref_sols)
    idx = []
    x_dist_Δ = []
    time_Δ = []

    for i in 1:n_maps
        for j in 1:n_x0s
            sol = sols[(i, j)]
            ref_sol = ref_sols[(i, j)]

            if sol.mcp_success && ref_sol.mcp_success
                push!(idx, (i, j))
                push!(x_dist_Δ, sol.x_dist - ref_sol.x_dist)
                push!(time_Δ, sol.time - ref_sol.time)
            end
        end
    end

    n_samples = length(idx)
    x_dist_Δ_CI = 1.96 * std(x_dist_Δ) / sqrt(n_samples)
    time_Δ_CI = 1.96 * std(time_Δ) / sqrt(n_samples)

    (; idx, x_dist_Δ, time_Δ, x_dist_Δ_CI, time_Δ_CI)
end