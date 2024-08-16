function get_sps_cons_3d(z, T, ego_polys, obs_polys, n_xu, n_per_col)
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    n_per_ego = n_per_col * n_obs
    n_per_t = n_per_col * n_obs * n_ego
    sides_per_obs = length(obs_polys[1].b)
    sides_per_ego = length(ego_polys[1].b)

    cons_sps = Num[]
    l_sps = Float64[]
    u_sps = Float64[]

    for t in 1:T
        xt = @view(z[(t-1)*n_xu+1:(t-1)*n_xu+6])
        yt = @view(z[T*n_xu+(t-1)*n_per_t+1:T*n_xu+(t-1)*n_per_t+n_per_t])

        for (i, Pe) in enumerate(ego_polys)
            yti = @view(yt[(i-1)*n_per_ego+1:(i-1)*n_per_ego+n_per_ego])
            Ve = Pe.V
            Vex = shift_to_3D(Ve, xt)

            for (j, Po) in enumerate(obs_polys)
                abc = @view(yti[(j-1)*n_per_col+1:(j-1)*n_per_col+3])
                d = yti[(j-1)*n_per_col+4]
                Vo = Po.V

                for k in 1:sides_per_ego
                    push!(cons_sps, abc' * Vex[k] + d) # >= 0 
                    push!(l_sps, 0.0)
                    push!(u_sps, Inf)
                end

                for k in 1:sides_per_obs
                    push!(cons_sps, -abc' * Vo[k] - d) # >= 0 
                    push!(l_sps, 0.0)
                    push!(u_sps, Inf)
                end

                push!(cons_sps, abc' * abc - 1) # = 0 
                push!(l_sps, -Inf)
                push!(u_sps, Inf)

            end
        end
    end
    (cons_sps, l_sps, u_sps)
end


function get_built_funs()
    xt = Symbolics.@variables(xt[1:6])[1] |> Symbolics.scalarize # state vector
    pa = Symbolics.@variables(pa[1:4])[1] |> Symbolics.scalarize # parameters of the hyperplane, pa = [a, b, c, d]
    V = Symbolics.@variables(V[1:3])[1] |> Symbolics.scalarize # representation of a vertex

    Vt = shift_to_3D([V], xt)[1]
    sep_ego_con = pa[1:3]'*Vt + pa[4] # ax+by+cz+d>=0, only for ego
    sep_obs_con = pa[1:3]'*V + pa[4] # ax+by+cz+d>=0, only for obs
    norm_con = pa[1:3]'*pa[1:3] - 1 # a²+b²+c²-1=0
    grad_ego_sep = Symbolics.gradient(sep_ego_con, [xt; pa])
    grad_obs_sep = Symbolics.gradient(sep_obs_con, pa)
    grad_norm = Symbolics.gradient(norm_con, pa)
    Hessian_ego_sep = Symbolics.jacobian(grad_ego_sep, [xt; pa])
    Hessian_obs_sep = Symbolics.jacobian(grad_obs_sep, pa)
    Hessian_norm = Symbolics.jacobian(grad_norm, pa)
    
    get_grad_ego_sep = Symbolics.build_function(grad_ego_sep, xt, pa, V; expression=Val(false))[2]
    get_grad_obs_sep = Symbolics.build_function(grad_obs_sep, pa, V; expression=Val(false))[2]
    get_grad_norm = Symbolics.build_function(grad_norm, pa; expression=Val(false))[2]
    get_Hessian_ego_sep = Symbolics.build_function(Hessian_ego_sep, xt, pa, V; expression=Val(false))[2]
    # get_Hessian_obs_sep = Symbolics.build_function(Hessian_obs_sep, pa, V; expression=Val(false))[2] # always zeros(4,4)
    Hessian_obs_sep = zeros(4, 4)
    # get_Hessian_norm = Symbolics.build_function(Hessian_norm, pa; expression=Val(false))[2] # always dia(2,2,2,0)
    Hessian_norm = Diagonal([2, 2, 2, 0])

    # return get_grad_ego_sep, get_grad_obs_sep, get_grad_norm, get_Hessian_ego_sep, get_Hessian_obs_sep, get_Hessian_norm
    return get_grad_ego_sep, get_grad_obs_sep, get_grad_norm, get_Hessian_ego_sep, Hessian_obs_sep, Hessian_norm
end

function sep_ego_con_val(xt, pa, V)
    Vt = shift_to_3D([V], xt)[1]
    pa[1:3]'*Vt + pa[4]
end

function sep_obs_con_val(pa, V)
    pa[1:3]'*V + pa[4]
end

function norm_con_val(pa)
    pa[1:3]'*pa[1:3] - 1
end

function setup_sep_planes_3d(
    ego_polys,
    obs_polys;
    T=1,
    dt=0.2,
    R_cost=1e-3 * I(6), # penality for control variables
    Q_cost=1e-3 * I(3), # penality for distance
    p1_max=500.0,
    p2_max=500.0,
    p3_max=500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=1.0,
    u4_max=π / 4,
    u5_max=π / 4,
    u6_max=π / 4,
)
    n_x = 12
    n_u = 6
    n_xu = n_x + n_u
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    # TODO not applicable if polys have different number of constraints
    n_side_ego = length(ego_polys[1].b)
    n_side_obs = length(obs_polys[1].b)
    # n_side_ego = maximum([length(i.b) for i in ego_polys]) # just for the buffer
    # n_side_obs = maximum([length(i.b) for i in obs_polys]) # just for the buffer
    n_per_col = 4 # for every collision, we need four parameters, a, b, c, and d, to represent a hyperplane ax + by + cz + d = 0 
    n_per_ego = n_per_col * n_obs
    n_per_t = n_per_col * n_obs * n_ego
    n_sps = (n_side_ego + n_side_obs + 1) * T # number of sep plane parameters

    n_dyn_cons = T * n_x
    n_env_cons = T * n_xu
    n_sep_cons = n_side_ego + n_side_obs + 1 # number of sep plane constraints
    xt_s2i, pa_s2i, dyn_cons_s2i, env_cons_s2i, sep_cons_s2i = get_sub2idxs((n_xu, T), (n_per_col, n_obs, n_ego, T), (n_dyn_cons), (n_env_cons), (n_sep_cons, n_obs, n_ego, T))

    z = Symbolics.@variables(z[1:n_xu*T+n_per_t*T])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:n_x])[1] |> Symbolics.scalarize

    # z consists of xt vectors and parameters of hyperplane, but the latter is not used in functions below
    cost_nom = f_3d(z, T, R_cost, Q_cost)
    cons_dyn = g_dyn_3d(z, x0, T, dt)
    cons_env = g_env_3d(z, T, p1_max, p2_max, p3_max, u1_max, u2_max, u3_max, u4_max, u5_max, u6_max)

    # check indexing consistency
    @assert length(dyn_cons_s2i) == length(cons_dyn)
    @assert length(env_cons_s2i) == length(cons_env)
    #Ve = ego_polys[1].V
    #Vos = map(obs_polys) do P
    #    hcat(P.V...)' |> collect
    #end

    #cons_sps2, l_sps2, u_sps2 = g_col_sps(z, T, Vos, Ve)
    cons_sps, l_sps, u_sps = get_sps_cons_3d(z, T, ego_polys, obs_polys, n_xu, n_per_col)

    cons_nom = [cons_dyn; cons_env] # cons_sps]

    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize
    λ_sps = Symbolics.@variables(λ_sps[1:n_sps])[1] |> Symbolics.scalarize

    θ = [z; λ_nom; λ_sps]

    lag = cost_nom - cons_nom' * λ_nom #- λ_sps' * cons_sps (to be filled later)
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; cons_nom; zeros(Num, n_sps)]
    F_nom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false))[2]

    n = length(F_nom)
    l = [fill(-Inf, length(grad_lag)); fill(-Inf, length(cons_dyn)); zeros(length(cons_env)); l_sps]
    u = [fill(+Inf, length(grad_lag)); fill(Inf, length(cons_dyn)); fill(Inf, length(cons_env)); u_sps]

    J_nom = Symbolics.sparsejacobian(F_nom, θ)
    (Jnom_rows, Jnom_cols, Jnom_vals) = findnz(J_nom)
    get_Jnom_vals! = Symbolics.build_function(Jnom_vals, z, x0, λ_nom; expression=Val(false))[2]

    xt = Symbolics.@variables(xt[1:6])[1] |> Symbolics.scalarize # state vector
    pa = Symbolics.@variables(pa[1:4])[1] |> Symbolics.scalarize # parameters of the hyperplane, pa = [a, b, c, d]
    get_grad_ego_sep, get_grad_obs_sep, get_grad_norm, get_Hessian_ego_sep, Hessian_obs_sep, Hessian_norm = get_built_funs()

    grad_xtpa_buf = zeros(10)
    grad_pa_buf = zeros(4)
    Hessian_xtpa_buf = zeros(10, 10)
    # Hessian_pa_buf = zeros(4, 4)

    function fill_F!(F, θ, x0)
        F .= 0.0
        @inbounds z = @view θ[[xt_s2i[:]; pa_s2i[:]]]
        @inbounds λ_nom = @view θ[[dyn_cons_s2i[:]; env_cons_s2i[:]]]
        F_nom!(F, z, x0, λ_nom)
        
        for t in 1:T
            xt_ind = xt_s2i[1:6, t] # only care about xt[1:3] and xt[4:6], i.e., position and orientation
            @inbounds xt = @view(θ[xt_ind])

            for (i, Pe) in enumerate(ego_polys)
                Ve = Pe.V
                m1 = length(Pe.b)

                for (j, Po) in enumerate(obs_polys)
                    pa_ind = pa_s2i[:, j, i, t]
                    @inbounds pa = @view θ[pa_ind]
                    Vo = Po.V
                    m2 = length(Po.b)
                    λ_ind = sep_cons_s2i[:, j, i, t]
                    @inbounds λ_sps = @view θ[λ_ind]
                    # @assert length(λ_sps) == m1+m2+1

                    for (k, V) in enumerate(Ve)
                        get_grad_ego_sep(grad_xtpa_buf, xt, pa, V)
                        F[[xt_ind; pa_ind]] .+= -λ_sps[k] * grad_xtpa_buf
                        F[λ_ind[k]] += sep_ego_con_val(xt, pa, V)
                    end
                    for (k, V) in enumerate(Vo)
                        get_grad_obs_sep(grad_pa_buf, pa, V) # notice there is no xt
                        F[pa_ind] .+= λ_sps[m1+k] * grad_pa_buf # notice there is no minus
                        F[λ_ind[m1+k]] += -sep_obs_con_val(pa, V)
                    end
                    get_grad_norm(grad_pa_buf, pa)
                    F[pa_ind] .+= -λ_sps[m1+m2+1] * grad_pa_buf
                    F[λ_ind[m1+m2+1]] += norm_con_val(pa)
                end
            end
        end
        

        nothing
    end

    Jnom_buf = zeros(length(Jnom_vals))
    function fill_J_vals!(J_vals, θ, x0)
        J_vals.nzval .= 1e-16 # clear
        @inbounds z = @view θ[[xt_s2i[:]; pa_s2i[:]]]
        @inbounds λ_nom = @view θ[[dyn_cons_s2i[:]; env_cons_s2i[:]]]
        get_Jnom_vals!(Jnom_buf, z, x0, λ_nom)
        J_vals .+= sparse(Jnom_rows, Jnom_cols, Jnom_buf, n, n)
        
        
        for t in 1:T
            xt_ind = xt_s2i[1:6, t] # only care about xt[1:3] and xt[4:6], i.e., position and orientation
            @inbounds xt = @view(θ[xt_ind])

            for (i, Pe) in enumerate(ego_polys)
                Ve = Pe.V
                m1 = length(Pe.b)

                for (j, Po) in enumerate(obs_polys)
                    pa_ind = pa_s2i[:, j, i, t]
                    @inbounds pa = @view θ[pa_ind]
                    Vo = Po.V
                    m2 = length(Po.b)
                    λ_ind = sep_cons_s2i[:, j, i, t]
                    @inbounds λ_sps = @view θ[λ_ind]
                    # @assert length(λ_sps) == m1+m2+1

                    for (k, V) in enumerate(Ve)
                        get_grad_ego_sep(grad_xtpa_buf, xt, pa, V)
                        get_Hessian_ego_sep(Hessian_xtpa_buf, xt, pa, V)

                        J_vals[[xt_ind;pa_ind], [xt_ind;pa_ind]] .+= -λ_sps[k] * Hessian_xtpa_buf
                        J_vals[[xt_ind;pa_ind], λ_ind[k]] .+= -grad_xtpa_buf
                        J_vals[λ_ind[k], [xt_ind;pa_ind]] .+= grad_xtpa_buf
                    end
                    for (k, V) in enumerate(Vo)
                        get_grad_obs_sep(grad_pa_buf, pa, V) # notice there is no xt

                        J_vals[pa_ind, pa_ind] .+= λ_sps[m1+k] * Hessian_obs_sep # do nothing
                        J_vals[pa_ind, λ_ind[m1+k]] .+= grad_pa_buf # notice there is no minus
                        J_vals[λ_ind[m1+k], pa_ind] .+= -grad_pa_buf # notice there is a minus
                    end
                    
                    get_grad_norm(grad_pa_buf, pa)

                    J_vals[pa_ind, pa_ind] .+= -λ_sps[m1+m2+1] * Hessian_norm # do nothing
                    J_vals[pa_ind, λ_ind[m1+m2+1]] .+= -grad_pa_buf
                    J_vals[λ_ind[m1+m2+1], pa_ind] .+= grad_pa_buf

                end
            end
        end

        nothing
    end

    J_example = sparse(Jnom_rows, Jnom_cols, ones(length(Jnom_cols)), n, n)
    for t in 1:T
        xt_ind = xt_s2i[1:n_x, t]

        for (i, Pe) in enumerate(ego_polys)
            for (j, Po) in enumerate(obs_polys)
                pa_ind = pa_s2i[:, j, i, t]
                sep_cons_ind = sep_cons_s2i[:, j, i, t]
                @inbounds J_example[[xt_ind;pa_ind;sep_cons_ind], [xt_ind;pa_ind;sep_cons_ind]] .+= 1.0
            end
        end
    end

    param = (;
        ego_polys,
        obs_polys,
        T,
        dt,
        R_cost,
        Q_cost,
        p1_max,
        p2_max,
        p3_max,
        u1_max,
        u2_max,
        u3_max,
        u4_max,
        u5_max,
        u6_max,
        xt_s2i,
        pa_s2i,
        dyn_cons_s2i,
        env_cons_s2i,
        sep_cons_s2i,
        n_per_col,
        n_per_ego,
        n_per_t,
        n_z=length(z)
    )

    return (; fill_F!,
        fill_J_vals!,
        J_example,
        l,
        u,
        param
    )
end

function visualize_sep_planes_3d(x0, T, ego_polys, obs_polys; n_per_col=4, fig=Figure(), ax3=LScene(fig[1, 1], scenekw=(camera=cam3d!, show_axis=true)), θ=[], is_displaying=true, is_showing_sep_plane=true)
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    n_per_t = n_per_col * n_obs * n_ego
    n_xu = 18
    xxts = Dict()
    abcdts = Dict()


    for i in 1:n_ego
        xx = x0[1:6]
        Aeb = shift_to_3D(ego_polys[i].A, ego_polys[i].b, xx)
        self_poly = ConvexPolygon3D(Aeb[1], Aeb[2])
        #plot!(ax, self_poly; color=:blue)
        plot_3D!(ax3, self_poly; color=:blue)

        for t in 1:T-1#5:1:T-1
            xxts[i, t] = Observable(x0[1:6])
            Aeb = @lift(shift_to_3D(ego_polys[i].A, ego_polys[i].b, $(xxts[i, t])))
            self_poly = @lift(ConvexPolygon3D($(Aeb)[1], $(Aeb)[2]))
            #plot!(ax, self_poly; color=:blue, linestyle=:dash)
            plot_3D!(ax3, self_poly; color=:blue, linestyle=:dash)
        end
        t = T
        xxts[i, t] = Observable(x0[1:6])
        Aeb = @lift(shift_to_3D(ego_polys[i].A, ego_polys[i].b, $(xxts[i, t])))
        self_poly = @lift(ConvexPolygon3D($(Aeb)[1], $(Aeb)[2]))
        #plot!(ax, self_poly; color=:blue, linewidth=3)
        plot_3D!(ax3, self_poly; color=:blue, linewidth=3)

        if is_showing_sep_plane
            for e in 1:n_obs
                abcdts[i, e] = Observable([1.0, 1.0, 1.0, -1.0])
                # ax + by + +cz + d = 0
                a = @lift($(abcdts[i, e])[1])
                b = @lift($(abcdts[i, e])[2])
                c = @lift($(abcdts[i, e])[3])
                d = @lift($(abcdts[i, e])[4])

                xs = LinRange(-10, 10, 3)
                ys = LinRange(-10, 10, 3)
                zs = @lift([(-($d)-($a)*x-($b)*y)/($c) for x in xs, y in ys])

                surface!(ax3, xs, ys, zs)
            end
        end
    end

    colors = [:red for _ in 1:n_obs]
    for (P, c) in zip(obs_polys, colors)
        plot_3D!(ax3, P; color=c)
        #plot!(ax, P; color=c)
    end

    function update_fig(θ)
        for i in 1:n_ego
            for t in 1:T #5:5:T
                xxts[i, t][] = copy(θ[(t-1)*18+1:(t-1)*18+6])
            end
            if is_showing_sep_plane
                for k in 1:n_obs
                    abcdts[i, k][] = copy(θ[n_xu*T+(T-1)*n_per_t+(k-1)*3+1:n_xu*T+(T-1)*n_per_t+(k-1)*3+4])
                end
            end
        end
    end

    if !isempty(θ)
        update_fig(θ)
    end

    if is_displaying
        display(fig)
    end

    (fig, update_fig, ax3)
end

function solve_prob_sep_planes_3d(prob, x0; θ0=nothing, is_displaying=true, sleep_duration=0.0)
    # (; F_both!, J_both, l, u, T, n_z, n_nom, n_sps, ego_polys, p1_max, p2_max, p3_max, u1_max, u2_max, u3_max, u4_max, u5_max, u6_max, n_xu, obs_polys, n_per_col) = prob
    (; fill_F!, fill_J_vals!, J_example, l, u, param) = prob
    param = prob.param
    T = param.T
    ego_polys = param.ego_polys
    obs_polys = param.obs_polys
    n_per_ego = param.n_per_ego
    n_per_col = param.n_per_col
    n_per_t = param.n_per_t
    n_z = param.n_z
    n_x = 12
    n_u = 6
    n_xu = n_x + n_u
    n = length(l)

    if is_displaying
        (fig, update_fig) = visualize_sep_planes_3d(x0, T, ego_polys, obs_polys; n_per_col)
    end

    if isnothing(θ0)
        θ0 = zeros(n)
        p = x0[1:3]
        mrp = x0[4:6]
        R = R_from_mrp(mrp)
        for t in 1:T
            θ0[(t-1)*n_xu+1:(t-1)*n_xu+12] = x0
            yt = @view(θ0[T*n_xu+(t-1)*n_per_t+1:T*n_xu+(t-1)*n_per_t+n_per_t])
            for (i, Pe) in enumerate(ego_polys)
                yti = @view(yt[(i-1)*n_per_ego+1:(i-1)*n_per_ego+n_per_ego])
                
                for (k, Po) in enumerate(obs_polys)
                    ce = R * Pe.c + p
                    co_to_ce = ce - Po.c
                    abc = co_to_ce / norm(co_to_ce)
                    d = -abc' * (ce+Po.c)/2
                    yti[(k-1)*n_per_col+1:(k-1)*n_per_col+3] .= abc
                    yti[(k-1)*n_per_col+4] = d
                end
            end
        end
        #for t in 1:T
        #    for i in 1:n_ego
        #        xxts[i, t][] = copy(θ0[(t-1)*n_xu+1:(t-1)*n_xu+6])
        #        for k in 1:n_obs
        #            abcdts[i, t, k][] = copy(θ0[n_xu*T+(t-1)*n_per_t+(k-1)*3+1:n_xu*T+(t-1)*n_per_t+(k-1)*3+3])
        #        end
        #    end
        #end
    end



    function F(n, θ, result)
        result .= 0.0
        fill_F!(result, θ, x0)
        if is_displaying
            update_fig(θ)
            if sleep_duration > 0
                sleep(sleep_duration)
            end
        end
        Cint(0)
    end

    J_vals = J_example
    J_col = J_vals.colptr[1:end-1]
    J_len = diff(J_vals.colptr)
    J_row = J_vals.rowval
    nnz_total = length(J_vals.nzval)

    function J(n, nnz, θ, col, len, row, data)
        @assert nnz == nnz_total
        data .= 0.0
        fill_J_vals!(J_vals, θ, x0)
        col .= J_col
        len .= J_len
        row .= J_row
        data .= J_vals.nzval
        Cint(0)
    end

    # force compilation
    buf = zeros(n)
    Jbuf = zeros(nnz_total)
    w = randn(length(θ0))
    F(n, w, buf)
    t=time()
    J(n, nnz_total, w, zero(J_col), zero(J_len), zero(J_row), Jbuf)
    println("update J takes ", time()-t, " seconds")

    # # check Jacobian
    # buf2 = zeros(n)
    # Jrows, Jcols, _ = findnz(prob.J_example)
    # Jnum = sparse(Jrows, Jcols, Jbuf)
    # Jnum2 = spzeros(n, n)
    # @info "Testing Jacobian accuracy numerically"
    # @showprogress for ni in 1:n
    #    wi = copy(w)
    #    wi[ni] += 1e-8
    #    F(n, wi, buf2)
    #    Jnum2[:, ni] = sparse((buf2 - buf) ./ 1e-8)
    #    if norm(buf2 - buf)>1e-5
    #     @infiltrate
    #    end
    # end
    # @info "Jacobian error is $(norm(Jnum2-Jnum))"

    PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")
    status, θ, info = PATHSolver.solve_mcp(
        F,
        J,
        l,
        u,
        θ0;
        silent=true,
        nnz=nnz_total,
        jacobian_structure_constant=true,
        output_linear_model="no",
        preprocess=1,
        output_warnings="no",
        jacobian_data_contiguous=true,
        cumulative_iteration_limit=100_000,
        convergence_tolerance=5e-4
    )

    fres = zeros(n)

    F(n, θ, fres)

    if is_displaying
        display(fig)
    end

    @inbounds z = @view(θ[1:n_z])

    (; status, info, θ, z, fres)
end
