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
            cost = f(res.z, param.T, param.Rf, param.Qf)
            final_pos = res.z[(param.T-1)*param.n_xu+1:(param.T-1)*param.n_xu+2]
            sols[i, j] = (; mcp_success, time, cost, final_pos, res)
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
            cost = f(res.z, param.T, param.Rf, param.Qf)
            final_pos = res.z[(param.T-1)*param.n_xu+1:(param.T-1)*param.n_xu+2]
            sols[i, j] = (; mcp_success, time, cost, final_pos, res)
            next!(p)
        end
    end
    sols
end


function multi_solve_dcol(ego_poly, x0s, maps, param)
    @argcheck length(x0s) == param.n_x0s
    @argcheck length(maps) == param.n_maps
    for map in maps
        @argcheck length(map) == param.n_obs
    end

    dcol_prob = PolyPlanning.setup_dcol(
        ego_poly;
        param.T,
        param.dt,
        param.Rf,
        param.Qf,
        param.u1_max,
        param.u2_max,
        param.u3_max,
        param.n_obs
    )

    sols = Dict()
    p = Progress(param.n_maps * param.n_x0s, dt=1.0)
    for (i, map) in enumerate(maps)
        for (j, x0) in enumerate(x0s)
            res = PolyPlanning.solve_dcol(dcol_prob, x0, map; is_displaying=false)
            mcp_success = res.status == PATHSolver.MCP_Solved
            time = res.info.total_time
            cost = f(res.z, param.T, param.Rf, param.Qf)
            final_pos = res.z[(param.T-1)*param.n_xu+1:(param.T-1)*param.n_xu+2]
            sols[i, j] = (; mcp_success, time, cost, final_pos, res)
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
            cost = f(res.z, param.T, param.Rf, param.Qf)
            final_pos = res.z[(param.T-1)*param.n_xu+1:(param.T-1)*param.n_xu+2]
            sols[i, j] = (; mcp_success, time, cost, final_pos, res)
            next!(p)
        end
    end
    sols
end

function load_all(exp_name, res_file_date, exp_file_date; is_loading_sep=false, is_loading_dcol=false, is_loading_kkt=false, data_dir="data")
    @info "Loading $exp_name exp results from $res_file_date for data from $exp_file_date..."
    our_file = jldopen("$data_dir/$(exp_name)_our_sols_$(res_file_date)_exp_$(exp_file_date).jld2", "r")
    our_sols = our_file["our_sols"]
    sep_sols = []
    kkt_sols = []
    if is_loading_sep
        sep_file = jldopen("$data_dir/$(exp_name)_sep_sols_$(res_file_date)_exp_$(exp_file_date).jld2", "r")
        sep_sols = sep_file["sep_sols"]
    end

    if is_loading_dcol
        dcol_file = jldopen("$data_dir/$(exp_name)_dcol_sols_$(res_file_date)_exp_$(exp_file_date).jld2", "r")
        dcol_sols = dcol_file["dcol_sols"]
    end

    if is_loading_kkt
        kkt_file = jldopen("$data_dir/$(exp_name)_kkt_sols_$(res_file_date)_exp_$(exp_file_date).jld2", "r")
        kkt_sols = kkt_file["kkt_sols"]
    end

    (our_sols, sep_sols, dcol_sols, kkt_sols)
end

function compute_all(ego_poly, x0s, maps, param; is_saving=false, exp_name="", is_running_sep=false, is_running_kkt=false, is_running_dcol=false, data_dir="data", date_now="", exp_file_date="")
    #n_maps = length(maps)
    #n_x0s = length(x0s)
    @info "Computing nonsmooth solutions..."
    start_t = time()
    our_sols = multi_solve_ours(ego_poly, x0s, maps, param)
    @info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."

    if is_saving
        jldsave("$data_dir/$(exp_name)_our_sols_$(date_now)_exp_$exp_file_date.jld2"; our_sols)
    end

    sep_sols = []
    if is_running_sep
        @info "Computing separating hyperplane solutions..."
        start_t = time()
        sep_sols = multi_solve_sep(ego_poly, x0s, maps, param)
        @info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."

        if is_saving
            jldsave("$data_dir/$(exp_name)_sep_sols_$(date_now)_exp_$exp_file_date.jld2"; sep_sols)
        end
    end

    dcol_sols = []
    if is_running_dcol
        @info "Computing DifferentiableCollisions.jl solutions..."
        start_t = time()
        dcol_sols = multi_solve_dcol(ego_poly, x0s, maps, param)
        @info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."

        if is_saving
            jldsave("$data_dir/$(exp_name)_dcol_sols_$(date_now)_exp_$exp_file_date.jld2"; dcol_sols)
        end
    end

    kkt_sols = []
    if is_running_kkt
        @info "Computing direct KKT solutions..."
        start_t = time()
        kkt_sols = multi_solve_kkt(ego_poly, x0s, maps, param)
        @info "Done! $(round(time() - start_t; sigdigits=3)) seconds elapsed."

        if is_saving
            jldsave("$data_dir/$(exp_name)_kkt_sols_$(date_now)_exp_$exp_file_date.jld2"; kkt_sols)
        end
    end

    (our_sols, sep_sols, dcol_sols, kkt_sols)
end

function load_experiment(name, date; data_dir="data")
    exp_file = jldopen("$data_dir/$(name)_exp_$date.jld2", "r")

    ego_poly = exp_file["ego_poly"]
    x0s = exp_file["x0s"]
    maps = exp_file["maps"]
    param = exp_file["param"]

    (; ego_poly, x0s, maps, param)
end


function process_into_bins(sols; global_success_radius=0.5)
    success = (idx=[], time=[], cost=[])
    #global_success = (idx=[], time=[], cost=[])
    fail = (idx=[], time=[], cost=[])

    for (i, sol) in sols
        if sol.mcp_success
            push!(success.idx, i)
            push!(success.time, sol.time)
            push!(success.cost, sol.cost)

            # for now, assume global success depends on local success
            #if sol.final_pos'sol.final_pos <= global_success_radius^2
            #    push!(global_success.idx, i)
            #    push!(global_success.time, sol.time)
            #    push!(global_success.cost, sol.cost)
            #end
        else
            push!(fail.idx, i)
            push!(fail.time, sol.time)
            push!(fail.cost, sol.cost)
        end
    end
    (; success, fail)
end

function find_bin_common(bin, ref_bin)
    idx = []
    time = []
    cost = []
    ref_time = []
    ref_cost = []


    for (i, id) in enumerate(bin.idx)
        i_ref = findfirst(x -> x == id, ref_bin.idx)

        if i_ref !== nothing
            push!(idx, id)
            push!(time, bin.time[i])
            push!(cost, bin.cost[i])
            push!(ref_time, ref_bin.time[i_ref])
            push!(ref_cost, ref_bin.cost[i_ref])
        end
    end

    (; idx, time, cost, ref_time, ref_cost)
end

function get_mean_std_CI(v)
    try
        m = mean(v)
        s = std(v)
        CI = 1.96 * std(v) / sqrt(length(v))
        return (mean=m, CI, std=s)
    catch e
        @error "empty v"
    end
end

function print_stats(bin, ref_bin, n_maps, n_x0s; name="bin", ref_name="ref_bin")

    common_success = PolyPlanning.find_bin_common(bin.success, ref_bin.success)
    common_fail = PolyPlanning.find_bin_common(bin.fail, ref_bin.fail)

    bin_success_percent = length(bin.success.idx) / (n_maps * n_x0s) * 100
    ref_success_percent = length(ref_bin.success.idx) / (n_maps * n_x0s) * 100
    both_success_percent = length(common_success.idx) / (n_maps * n_x0s) * 100
    bin_fail_percent = length(bin.fail.idx) / (n_maps * n_x0s) * 100
    ref_fail_percent = length(ref_bin.fail.idx) / (n_maps * n_x0s) * 100
    both_fail_percent = length(common_fail.idx) / (n_maps * n_x0s) * 100

    println("           $name     $ref_name     both    both(ct)")
    println("success    $(round(bin_success_percent; sigdigits=3))%    $(round(ref_success_percent; sigdigits=3))%    $(round(both_success_percent; sigdigits=3))%    $(length(common_success.idx))")
    println("fail       $(round(bin_fail_percent; sigdigits=3))%    $(round(ref_fail_percent; sigdigits=3))%    $(round(both_fail_percent; sigdigits=3))%    $(length(common_fail.idx))")
    if length(common_success.idx) > 0
        println("$name vs $ref_name successes compared:")
        print_table(common_success.time, common_success.ref_time, common_success.cost, common_success.ref_cost; name, ref_name)
    end
    #if length(common_fail.idx) > 0
    #println("fails compared:")
    #print_table(common_fail.time, common_fail.ref_time, common_fail.cost, common_fail.ref_cost; name, ref_name)
    #end
end

function print_table(time, ref_time, cost, ref_cost; name="bin", ref_name="ref_bin")

    abs_Δ_time = get_mean_std_CI(time - ref_time)
    rel_Δ_time = abs_Δ_time.mean ./ mean(ref_time) * 100
    rel_Δ_time_CI = abs_Δ_time.CI ./ mean(ref_time) * 100
    abs_Δ_cost = get_mean_std_CI(cost - ref_cost)
    rel_Δ_cost = abs_Δ_cost.mean ./ mean(ref_cost) * 100
    rel_Δ_cost_CI = abs_Δ_cost.CI ./ mean(ref_cost) * 100

    println("              mean    CI     std")
    println("time abs Δ    $(round.(abs_Δ_time.mean; sigdigits=2)) ± $(round.(abs_Δ_time.CI; sigdigits=2))  $(round.(abs_Δ_time.std; sigdigits=2))")
    println("time rel Δ    $(round.(rel_Δ_time; sigdigits=3))% ± $(round.(rel_Δ_time_CI; sigdigits=3))%  ")
    println("cost abs Δ    $(round.(abs_Δ_cost.mean; sigdigits=3)) ± $(round.(abs_Δ_cost.CI; sigdigits=3))  $(round.(abs_Δ_cost.std; sigdigits=3))")
    println("cost rel Δ    $(round.(rel_Δ_cost; sigdigits=3))% ± $(round.(rel_Δ_cost_CI; sigdigits=3))%  ")

    time_stats = get_mean_std_CI(time)
    ref_time_stats = get_mean_std_CI(ref_time)
    cost_stats = get_mean_std_CI(cost)
    ref_cost_stats = get_mean_std_CI(ref_cost)
    println("mean   $name vs $ref_name")
    println("time   $(round.(time_stats.mean; sigdigits=2))    $(round.(ref_time_stats.mean; sigdigits=2))  ")
    println("cost   $(round.(cost_stats.mean; sigdigits=2))    $(round.(ref_cost_stats.mean; sigdigits=2))  ")
end
#function print_table(bin, n_maps, n_x0s; title="our")
#our_success_ratio = length(bin.successes.idx) / (n_maps * n_x0s)
#    #our_global_success_ratio = length(bin.global_success.idx) / (n_maps * n_x0s)
#    our_fail_ratio = length(bin.fails.idx) / (n_maps * n_x0s)

#    println("             local       fails")
#    println("$title ratio    $(round(our_success_ratio*100; sigdigits=2))%      $(round(our_fail_ratio*100; sigdigits=2))%")
#    println("          (mean, CI)  (mean, CI)")
#    println("$title time $(round.(get_mean_CI(bin.successes.time); sigdigits=2)) $(round.(get_mean_CI(bin.fails.time); sigdigits=2))")
#    println("$title cost $(round.(get_mean_CI(bin.successes.cost); sigdigits=2)) $(round.(get_mean_CI(bin.fails.cost); sigdigits=2))")
#end


function visualize_multi(x0s, maps, sols, bins, T, ego_poly; n_rows=1, n_cols=1, is_displaying=true, type="nonsmooth")

    for k in 1:n_rows*n_cols:length(bins.idx)
        counter = 0
        fig = Figure()
        idx_range = k:min(k + n_rows * n_cols - 1, length(bins.idx))

        for (idx, time, cost) in zip(bins.idx[idx_range], bins.time[idx_range], bins.cost[idx_range])
            sols
            map_idx = idx[1]
            x0_idx = idx[2]
            x0 = x0s[x0_idx]
            map = maps[map_idx]
            sol = sols[(map_idx, x0_idx)]
            counter += 1

            # ind2sub...
            i = Int(floor((counter - 1) / n_cols)) + 1
            j = Int((counter - 1) % n_cols) + 1
            #@info "$i, $j"
            ax = Axis(fig[i, j], aspect=DataAspect())
            ax.title = "$type\nmaps[$(map_idx)], x0s[$(x0_idx)] = $(round.(x0[1:3];sigdigits=2)), $(sol.mcp_success ? "success" : "FAIL")\ntime = $(round(sol.time; sigdigits=2)) s, cost = $(round(sol.cost; sigdigits=2)), pT = $(round.(sol.final_pos; sigdigits=2)) "

            if type == "nonsmooth"
                if sol.mcp_success
                    (fig, update_fig, ax) = visualize_quick(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                else
                    (fig, update_fig, ax) = visualize_quick(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                    lines!(ax, [-5, 5], [-5, 5]; color=:red, linewidth=10)
                end
            elseif type == "sep_planes"
                if sol.mcp_success
                    (fig, update_fig, ax) = visualize_sep_planes(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                else
                    (fig, update_fig, ax) = visualize_sep_planes(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                    lines!(ax, [-5, 5], [-5, 5]; color=:red, linewidth=10)
                end
            elseif type == "dcol"
                if sol.mcp_success
                    (fig, update_fig, ax) = visualize_dcol(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                else
                    (fig, update_fig, ax) = visualize_dcol(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                    lines!(ax, [-5, 5], [-5, 5]; color=:red, linewidth=10)
                end
            elseif type == "direct_kkt"
                if sol.mcp_success
                    (fig, update_fig, ax) = visualize_direct_kkt(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                else
                    (fig, update_fig, ax) = visualize_direct_kkt(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                    lines!(ax, [-5, 5], [-5, 5]; color=:red, linewidth=10)
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


function visualize_multi(x0s, maps, sols, T, ego_poly; n_rows=1, n_cols=1, is_displaying=true, type="nonsmooth")
    n_maps = length(maps)
    n_x0s = length(x0s)
    @argcheck n_maps * n_x0s == length(sols)

    n_rows = min(n_maps, n_rows)
    n_cols = min(n_x0s, n_cols)

    for maps_idx_begin in 1:n_rows:n_maps
        for x0_idx_begin in 1:n_cols:n_x0s
            fig = Figure()
            #empty!(fig.scene) # clear figure, interactivity doesn't after first page with this

            for i in 1:n_rows
                for j in 1:n_cols
                    map_idx = maps_idx_begin - 1 + i
                    x0_idx = x0_idx_begin - 1 + j

                    if map_idx <= n_maps && x0_idx <= n_x0s
                        x0 = x0s[x0_idx]
                        map = maps[map_idx]
                        sol = sols[(map_idx, x0_idx)]
                        ax = Axis(fig[i, j], aspect=DataAspect())

                        ax.title = "$type\nmaps[$(map_idx)], x0s[$(x0_idx)] = $(round.(x0[1:3];sigdigits=2)), $(sol.mcp_success ? "success" : "FAIL")\ntime = $(round(sol.time; sigdigits=2)) s, cost = $(round(sol.cost; sigdigits=2)), pT = $(round.(sol.final_pos; sigdigits=2)) "

                        if type == "nonsmooth"
                            if sol.mcp_success
                                (fig, update_fig, ax) = visualize_quick(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                            else
                                (fig, update_fig, ax) = visualize_quick(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                                lines!(ax, [-5, 5], [-5, 5]; color=:red, linewidth=10)
                            end
                        elseif type == "sep_planes"
                            if sol.mcp_success
                                (fig, update_fig, ax) = visualize_sep_planes(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                            else
                                (fig, update_fig, ax) = visualize_sep_planes(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                                lines!(ax, [-5, 5], [-5, 5]; color=:red, linewidth=10)
                            end
                        elseif type == "dcol"
                            if sol.mcp_success
                                (fig, update_fig, ax) = visualize_dcol(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                            else
                                (fig, update_fig, ax) = visualize_dcol(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                                lines!(ax, [-5, 5], [-5, 5]; color=:red, linewidth=10)
                            end
                        elseif type == "direct_kkt"
                            if sol.mcp_success
                                (fig, update_fig, ax) = visualize_direct_kkt(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                            else
                                (fig, update_fig, ax) = visualize_direct_kkt(x0, T, ego_poly, map; fig, ax, sol.res.θ, is_displaying=false)
                                lines!(ax, [-5, 5], [-5, 5]; color=:red, linewidth=10)
                            end
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