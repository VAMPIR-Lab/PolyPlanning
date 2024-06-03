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
        Q=0.0 * PolyPlanning.I(2), # disabled final cost
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
            final_pos = res.z[(param.T-1)*param.n_xu+1:(param.T-1)*param.n_xu+2]
            x_dist = final_pos[1]
            sols[i, j] = (; mcp_success, time, x_dist, final_pos, res)
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
            final_pos = res.z[(param.T-1)*param.n_xu+1:(param.T-1)*param.n_xu+2]
            x_dist = final_pos[1]
            sols[i, j] = (; mcp_success, time, x_dist, final_pos, res)
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
            final_pos = res.z[(param.T-1)*param.n_xu+1:(param.T-1)*param.n_xu+2]
            x_dist = final_pos[1]
            sols[i, j] = (; mcp_success, time, x_dist, final_pos, res)
            next!(p)
        end
    end
    sols
end

function load_all(exp_name, res_file_date, exp_file_date; is_loading_sep=false, is_loading_kkt=false, data_dir="data")
    @info "Loading $exp_name exp results from $res_file_date for data from $exp_file_date..."
    our_file = jldopen("$data_dir/$(exp_name)_our_sols_$(res_file_date)_exp_$(exp_file_date).jld2", "r")
    our_sols = our_file["our_sols"]
    sep_sols = []
    kkt_sols = []
    if is_loading_sep
        sep_file = jldopen("$data_dir/$(exp_name)_sep_sols_$(res_file_date)_exp_$(exp_file_date).jld2", "r")
        sep_sols = sep_file["sep_sols"]
    end

    if is_loading_kkt
        kkt_file = jldopen("$data_dir/$(name)_kkt_sols_$(res_date)_exp_$(exp_date).jld2", "r")
        kkt_sols = kkt_file["kkt_sols"]
    end

    (our_sols, sep_sols, kkt_sols)
end

function compute_all(ego_poly, x0s, maps, param; is_saving=false, exp_name="", is_running_sep=false, is_running_kkt=false, data_dir="data", date_now="", exp_file_date="")
    n_maps = length(maps)
    n_x0s = length(x0s)
    @info "Computing nonsmooth solutions..."
    start_t = time()
    our_sols = multi_solve_ours(ego_poly, x0s, maps, param)
    @info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
    our_only_mcp_success = filter_by_mcp_success(our_sols)
    @info "$(length(our_only_mcp_success.idx)/(n_maps*n_x0s)*100)% nonsmooth success rate"

    if is_saving
        jldsave("$data_dir/$(exp_name)_our_sols_$(date_now)_exp_$exp_file_date.jld2"; our_sols)
    end

    sep_sols = []
    kkt_sols = []

    if is_running_sep
        @info "Computing separating hyperplane solutions..."
        start_t = time()
        sep_sols = multi_solve_sep(ego_poly, x0s, maps, param)
        @info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
        sep_only_mcp_success = filter_by_mcp_success(sep_sols)
        @info "$(length(sep_only_mcp_success.idx)/(n_maps*n_x0s)*100)% sep mcp success"

        if is_saving
            jldsave("$data_dir/$(exp_name)_sep_sols_$(date_now)_exp_$exp_file_date.jld2"; sep_sols)
        end
    end

    if is_running_kkt
        @info "Computing direct KKT solutions..."
        start_t = time()
        kkt_sols = multi_solve_kkt(ego_poly, x0s, maps, param)
        @info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."
        kkt_only_mcp_success = filter_by_mcp_success(kkt_sols)
        @info "$(length(kkt_only_mcp_success.idx)/(n_maps*n_x0s)*100)% kkt mcp success rate"

        if is_saving
            jldsave("$data_dir/$(exp_name)_kkt_sols_$(date_now)_exp_$exp_file_date.jld2"; kkt_sols)
        end
    end

    (our_sols, sep_sols, kkt_sols)
end

function load_experiment(name, date; data_dir="data")
    exp_file = jldopen("$data_dir/$(name)_exp_$date.jld2", "r")

    ego_poly = exp_file["ego_poly"]
    x0s = exp_file["x0s"]
    maps = exp_file["maps"]
    param = exp_file["param"]

    (; ego_poly, x0s, maps, param)
end

function filter_by_mcp_success(sols)
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

function filter_by_task_success(sols; task_radius=0.5)
    idx = []
    times = []
    x_dists = []

    for (i, sol) in sols
        if sol.final_pos'sol.final_pos <= task_radius^2
            push!(idx, i)
            push!(times, sol.time)
            push!(x_dists, sol.x_dist)
        end
    end
    (; idx, times, x_dists)
end

function compute_sols_Δ_mcp(n_maps, n_x0s, sols, ref_sols)
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

function compute_sols_Δ_task(n_maps, n_x0s, sols, ref_sols; task_radius=0.5)
    idx = []
    x_dist_Δ = []
    time_Δ = []

    for i in 1:n_maps
        for j in 1:n_x0s
            sol = sols[(i, j)]
            ref_sol = ref_sols[(i, j)]

            if sol.final_pos'sol.final_pos <= task_radius^2 && ref_sol.final_pos'ref_sol.final_pos <= task_radius^2
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

function visualize_multi_quick(x0s, maps, sols, T, ego_poly; n_rows=1, n_cols=1, is_displaying=true, title_prefix="")
    n_maps = length(maps)
    n_x0s = length(x0s)
    @argcheck n_maps * n_x0s == length(sols)

    n_rows = min(n_maps, n_rows)
    n_cols = min(n_x0s, n_cols)

    fig = Figure()

    for maps_idx_begin in 1:n_rows:n_maps
        for x0_idx_begin in 1:n_cols:n_x0s
            empty!(fig.scene) # clear figure
            for i in 1:n_rows
                for j in 1:n_cols
                    maps_idx = maps_idx_begin - 1 + i
                    x0_idx = x0_idx_begin - 1 + j

                    if maps_idx <= n_maps && x0_idx <= n_x0s
                        x0 = x0s[x0_idx]
                        map = maps[maps_idx]
                        sol = sols[(maps_idx, x0_idx)]
                        ax = Axis(fig[i, j], aspect=DataAspect())

                        ax.title = "$title_prefix\nmaps[$(maps_idx)], x0s[$(x0_idx)] = $(round.(x0[1:3];sigdigits=2))\nmcp $(sol.mcp_success ? "success" : "FAIL"), $(round(sol.time; sigdigits=2)) s, pT = $(round.(sol.final_pos; sigdigits=2)) "

                        if sol.mcp_success
                            (fig, update_fig, ax) = visualize_quick(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                        else
                            (fig, update_fig, ax) = visualize_quick(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                            lines!(ax, [-5, 5], [-5, 5]; color=:red, linewidth=10)
                        end
                    end
                end
            end

            if is_displaying
                display(fig)
                @info "Enter nothing to continue, enter anything to stop..."
                if readline() != ""
                    return
                end
            end
        end
    end
end