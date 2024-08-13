# given a hyperplane a1*x + a2*y + b = 0, return line segment inside the limits
function get_line_segment_3d(a1, a2, b; lim_x=[-10,10], lim_y=[-10,10])
    if a1 == 0
        xs = lim_x
        ys = [-b/a2, -b/a2]
    else
        ys = lim_y
        xs = (-a2 * ys .-b) ./ a1
    end
    return xs, ys
end

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

        for (i, Pi) in enumerate(ego_polys)
            yti = @view(yt[(i-1)*n_per_ego+1:(i-1)*n_per_ego+n_per_ego])
            Ve = Pi.V
            Vex = shift_to_3D(Ve, xt)

            for (k, Pk) in enumerate(obs_polys)
                abc = @view(yti[1:3])
                d = yti[4]
                Vo = Pk.V

                for j in 1:sides_per_ego
                    push!(cons_sps, abc' * Vex[j] + d) # >= 0 
                    push!(l_sps, 0.0)
                    push!(u_sps, Inf)
                end

                for j in 1:sides_per_obs
                    push!(cons_sps, -abc' * Vo[j] - d) # >= 0 
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
    xdim = 12
    udim = 6
    n_xu = xdim + udim
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    # n_side_ego = length(ego_polys[1].b)
    # n_side_obs = length(obs_polys[1].b)
    # n_side_ego = maximum([length(i.b) for i in ego_polys]) # just for the buffer
    # n_side_obs = maximum([length(i.b) for i in obs_polys]) # just for the buffer
    n_per_col = 4 # for every collision, we need four parameters, a, b, c, and d, to represent a hyperplane ax + by + cz + d = 0 
    n_per_ego = n_per_col * n_obs
    n_per_t = n_per_col * n_obs * n_ego

    z = Symbolics.@variables(z[1:n_xu*T+n_per_t*T])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:xdim])[1] |> Symbolics.scalarize

    cost_nom = f_3d(z, T, R_cost, Q_cost)
    cons_dyn = g_dyn_3d(z, x0, T, dt)
    cons_env = g_env_3d(z, T, p1_max, p2_max, p3_max, u1_max, u2_max, u3_max, u4_max, u5_max, u6_max)

    #Ve = ego_polys[1].V
    #Vos = map(obs_polys) do P
    #    hcat(P.V...)' |> collect
    #end

    #cons_sps2, l_sps2, u_sps2 = g_col_sps(z, T, Vos, Ve)
    cons_sps, l_sps, u_sps = get_sps_cons_3d(z, T, ego_polys, obs_polys, n_xu, n_per_col)

    cons_nom = [cons_dyn; cons_env; cons_sps]

    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize

    θ = [z; λ_nom]

    lag = cost_nom - cons_nom' * λ_nom
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; cons_nom]
    F_nom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false), parallel=Symbolics.SerialForm())[2]

    l = [fill(-Inf, length(grad_lag)); fill(-Inf, length(cons_dyn)); zeros(length(cons_env)); l_sps]
    u = [fill(+Inf, length(grad_lag)); fill(Inf, length(cons_dyn)); fill(Inf, length(cons_env)); u_sps]

    J_nom = Symbolics.sparsejacobian(F_nom, θ)
    (J_rows_nom, J_cols_nom, J_vals) = findnz(J_nom)
    J_vals_nom! = Symbolics.build_function(J_vals, z, x0, λ_nom; expression=Val(false), parallel=Symbolics.SerialForm())[2]

    function F_both!(F, z_local, x0_local, λ_nom_local)
        F .= 0.0
        F_nom!(F, z_local, x0_local, λ_nom_local)
        nothing
    end

    function J_both_vals!(J_vals, z_local, x0_local, λ_nom_local)
        J_vals .= 0.0
        J_vals_nom!(J_vals, z_local, x0_local, λ_nom_local)
        nothing
    end

    return (; F_both!,
        J_both=(J_rows_nom, J_cols_nom, J_both_vals!),
        l,
        u,
        T,
        n_z=length(z),
        n_nom=length(λ_nom),
        ego_polys,
        p1_max,
        p2_max,
        p3_max,
        u1_max,
        u2_max,
        u3_max,
        u4_max,
        u5_max,
        u6_max,
        n_xu,
        obs_polys,
        n_per_col,
        n_per_ego,
        n_per_t
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
    (; F_both!, J_both, l, u, T, n_z, n_nom, ego_polys, p1_max, p2_max, p3_max, u1_max, u2_max, u3_max, u4_max, u5_max, u6_max, n_xu, obs_polys, n_per_col) = prob

    J_rows, J_cols, J_vals! = J_both
    nnz_total = length(J_rows)
    n = length(l)
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    n_per_ego = n_per_col * n_obs
    n_per_t = n_per_col * n_obs * n_ego

    @assert n == n_z + n_nom "did you forget to update l/u"

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
                    #Vo = hcat(Pk.V...)' |> collect
                    #Vo_center = sum(Vo; dims=1) ./ size(Vo, 1)
                    #a = x0[1:2] - Vo_center[1:2]
                    #z = (x0[1:2] + Vo_center[1:2]) / 2
                    #b = -a'z
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

    J_shape = sparse(J_rows, J_cols, Vector{Cdouble}(undef, nnz_total), n, n)
    J_col = J_shape.colptr[1:end-1]
    J_len = diff(J_shape.colptr)
    J_row = J_shape.rowval

    function F(n, θ, result)
        result .= 0.0
        @inbounds z = θ[1:n_z]
        @inbounds λ_nom = θ[n_z+1:n_z+n_nom]
        F_both!(result, z, x0, λ_nom)
        if is_displaying
            update_fig(θ)
            if sleep_duration > 0
                sleep(sleep_duration)
            end
        end
        Cint(0)
    end
    function J(n, nnz, θ, col, len, row, data)
        @assert nnz == nnz_total
        data .= 0.0
        @inbounds z = θ[1:n_z]
        @inbounds λ_nom = θ[n_z+1:n_z+n_nom]
        J_vals!(data, z, x0, λ_nom)
        col .= J_col
        len .= J_len
        row .= J_row
        Cint(0)
    end

    # force compilation
    buf = zeros(n)
    Jbuf = zeros(nnz_total)
    w = randn(length(θ0))
    F(n, w, buf)
    J(n, nnz_total, w, zero(J_col), zero(J_len), zero(J_row), Jbuf)

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
