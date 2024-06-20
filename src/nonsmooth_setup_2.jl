function gen_LP_data(xt, A_ego::AbstractArray{T}, b_ego, centr_ego, A_obs, b_obs, centr_obs) where {T}
    A = [A_ego A_ego*centr_ego+b_ego
        A_obs A_obs*centr_obs+b_obs]
    b = [b_ego; b_obs]
    q = [0, 0, 1.0]
    (A, b, q)
end

function is_ass_feasible(ass, m1, m2)
    if_ego = false
    if_obs = false
    for i in ass
        if i ∈ [i for i in 1:m1]
            if_ego = true
        end
        if i ∈ [i for i in m1+1:m1+m2]
            if_obs = true
        end
    end
    return length(ass) <= 3 && if_ego && if_obs
end

function g_col_single(xt, A_ego, b_ego, centr_ego, A_obs, b_obs, centr_obs)
    sds = Dict()
    Aex, bex = shift_to(A_ego, b_ego, xt)
    R = [cos(xt[3]) sin(xt[3])
        -sin(xt[3]) cos(xt[3])]
    centroidex = xt[1:2] + R * centr_ego
    AA, bb, qq = gen_LP_data(xt, Aex, bex, centroidex, A_obs, b_obs, centr_obs)
    m1 = length(bex)
    m2 = length(b_obs)
    all_active_inds = collect(1:m1+m2)
    Itr = powerset(all_active_inds) |> collect

    for active_inds in Itr
        if !is_ass_feasible(active_inds, m1, m2)
            continue
        end

        try
            AA_active = collect(AA[active_inds, :])
            bb_active = collect(bb[active_inds])
            # TODO what if not unique primal? Need to resolve
            if length(active_inds) == 3
                zz = -AA_active \ bb_active
            else
                # if linear system is underdetermined, use minimum norm solution (calculated by right inverse)
                # Note: every solution has the same sd, i.e., zz[3], but different zz[1:2]
                zz = -AA_active' * ((AA_active * AA_active') \ bb_active)
            end
            sd = zz[3]
            sds[active_inds] = sd
        catch err
            if err isa LinearAlgebra.SingularException
                continue
            else
                @warn(err)
            end
        end
    end
    sds
end

function setup_nonsmooth(
    ego_polys,
    obs_polys;
    T=2,
    dt=0.2,
    R_cost=1e-3 * I(3),
    Q_cost=1e-3 * I(2),
    p1_max=500.0,
    p2_min=-500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=π / 4
)

    # problem dimensions
    n_x = 6
    n_u = 3
    n_xu = n_x + n_u
    n_z = T * n_xu
    n_ego = length(ego_polys)
    n_obs = length(obs_polys)
    n_side_ego = length(ego_polys[1].b)
    n_side_obs = length(obs_polys[1].b)
    n_dyn_cons = T * n_x
    n_env_cons = T * n_xu
    combin_2_from_n = n::Int -> n * (n - 1) ÷ 2
    n_sd_cons_per_col = (n_side_ego * n_side_obs + combin_2_from_n(n_side_ego) * n_side_obs + n_side_ego * combin_2_from_n(n_side_obs))
    n_sd_cons = T * n_ego * n_obs * n_sd_cons_per_col

    # assume number of sides are the same for all ego polys, obs polys
    @assert all([length(ego_polys[i].b) for i in 1:n_ego] .== n_side_ego)
    @assert all([length(obs_polys[i].b) for i in 1:n_obs] .== n_side_obs)

    # θ indexing
    #z_s2i, dyn_cons_s2i, env_cons_s2i, sd_cons_s2i = get_sub2idxs((n_xu, T), (n_dyn_cons), (n_env_cons), (n_sd_cons_per_col, n_ego, n_obs, T))
    # hacky
    z_s2i, dyn_cons_s2i, env_cons_s2i, sd_cons_s2i = get_sub2idxs((n_xu, T), (n_dyn_cons), (n_env_cons), (2, n_ego, n_obs, T))

    z = Symbolics.@variables(z[1:n_z])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:n_x])[1] |> Symbolics.scalarize
    xt = Symbolics.@variables(xt[1:n_x])[1] |> Symbolics.scalarize

    cost = f(z, T, R_cost, Q_cost)
    dyn_cons = g_dyn(z, x0, T, dt)
    env_cons = g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)

    # indexing consistency
    @assert length(z_s2i) == n_z
    @assert length(dyn_cons_s2i) == length(dyn_cons)
    @assert length(env_cons_s2i) == length(env_cons)
    # disabled for hack
    #@assert length(sd_cons_s2i) == n_sd_cons

    sds = map(ego_polys) do Pe
        map(obs_polys) do Po
            A_ego = collect(Pe.A)
            A_obs = collect(Po.A)
            b_ego = Pe.b
            b_obs = Po.b
            V_ego = Pe.V
            V_obs = Po.V
            centr_ego = sum(V_ego) / length(V_ego)
            centr_obs = sum(V_obs) / length(V_obs)
            g_col_single(xt, A_ego, b_ego, centr_ego, A_obs, b_obs, centr_obs)
        end
    end

    # instead of something smarter:
    #@infiltrate
    sd_cons_hack =
        mapreduce(vcat, 1:T) do t
            [
                substitute(sds[1][1][[1, 2, 8]], xt .=> z[z_s2i[1:n_x, t]])
                substitute(sds[1][1][[1, 4, 8]], xt .=> z[z_s2i[1:n_x, t]])
            ]
        end

    #@infiltrate

    all_cons = [dyn_cons; env_cons; sd_cons_hack]
    λ_all = Symbolics.@variables(λ_nom[1:length(all_cons)])[1] |> Symbolics.scalarize

    θ = [z; λ_all]

    lag = cost - all_cons' * λ_all
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; all_cons]

    n = length(F_nom)
    l = zeros(n)
    u = zeros(n)

    # gradient = 0
    l[z_s2i[:]] .= -Inf
    u[z_s2i[:]] .= Inf

    # dynamics constraints
    l[dyn_cons_s2i[:]] .= -Inf
    u[dyn_cons_s2i[:]] .= Inf

    # environment constraints
    l[env_cons_s2i[:]] .= 0.0
    u[env_cons_s2i[:]] .= Inf

    # sd constraints
    l[sd_cons_s2i[:]] .= 0.0
    u[sd_cons_s2i[:]] .= Inf

    get_Fnom! = Symbolics.build_function(F_nom, z, x0, λ_all; expression=Val(false))[2]

    J_nom = Symbolics.sparsejacobian(F_nom, θ)
    (J_rows_nom, J_cols_nom, J_vals) = findnz(J_nom)
    J_vals_nom! = Symbolics.build_function(J_vals, z, x0, λ_all; expression=Val(false), parallel=Symbolics.SerialForm())[2]


    function F_both!(F, z_local, x0_local, λ_nom_local)
        F .= 0.0
        get_Fnom!(F, z_local, x0_local, λ_nom_local)
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
        ego_polys,
        obs_polys,
        p1_max,
        p2_min,
        z_s2i,
        dyn_cons_s2i,
        env_cons_s2i,
        sd_cons_s2i
    )
end


function visualize_nonsmooth(x0, T, ego_polys, obs_polys; fig=Figure(), ax=Axis(fig[1, 1], aspect=DataAspect()), θ=[], is_displaying=true, is_newsd=false)
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    xxts = Dict()

    for i in 1:n_ego
        xx = x0[1:3]
        Aeb = shift_to(ego_polys[i].A, ego_polys[i].b, xx)
        self_poly = ConvexPolygon2D(Aeb[1], Aeb[2])
        plot!(ax, self_poly; color=:blue)
        for t in 1:T-1#5:1:T-1
            xxts[i, t] = Observable(x0[1:3])
            Aeb = @lift(shift_to(ego_polys[i].A, ego_polys[i].b, $(xxts[i, t])))
            self_poly = @lift(ConvexPolygon2D($(Aeb)[1], $(Aeb)[2]))
            plot!(ax, self_poly; color=:blue, linestyle=:dash)
        end
        t = T
        xxts[i, t] = Observable(x0[1:3])
        Aeb = @lift(shift_to(ego_polys[i].A, ego_polys[i].b, $(xxts[i, t])))
        self_poly = @lift(ConvexPolygon2D($(Aeb)[1], $(Aeb)[2]))
        plot!(ax, self_poly; color=:blue, linewidth=3)
    end

    colors = [:red for _ in 1:n_obs]
    for (P, c) in zip(obs_polys, colors)
        plot!(ax, P; color=c)
    end

    function update_fig(θ)
        for i in 1:n_ego
            for t in 1:T #5:5:T
                xxts[i, t][] = copy(θ[(t-1)*9+1:(t-1)*9+6])
            end
        end
    end

    if !isempty(θ)
        update_fig(θ)
    end

    if is_displaying
        display(fig)
    end

    (fig, update_fig, ax)
end

function solve_nonsmooth(prob, x0; θ0=nothing, is_displaying=true, sleep_duration=0.0)
    (; F_both!, J_both, l, u, T, ego_polys, obs_polys, p1_max, p2_min, z_s2i, dyn_cons_s2i, env_cons_s2i, sd_cons_s2i) = prob

    n_x = 6
    n_u = 3
    n_xu = n_x + n_u
    J_rows, J_cols, J_vals! = J_both
    nnz_total = length(J_rows)
    n = length(l)
    #n_obs = length(obs_polys)
    #n_ego = length(ego_polys)
    #@assert n == n_z + n_nom "did you forget to update l/u"
    λ_all_s2i = [dyn_cons_s2i...; env_cons_s2i...; sd_cons_s2i...]

    if is_displaying
        (fig, update_fig) = visualize_nonsmooth(x0, T, ego_polys, obs_polys)
    end

    if isnothing(θ0)
        θ0 = zeros(n)
        for t in 1:T
            θ0[(t-1)*n_xu+1:(t-1)*n_xu+n_x] = x0
        end
    end

    J_shape = sparse(J_rows, J_cols, Vector{Cdouble}(undef, nnz_total), n, n)
    J_col = J_shape.colptr[1:end-1]
    J_len = diff(J_shape.colptr)
    J_row = J_shape.rowval

    function F(n, θ, result)
        result .= 0.0
        @inbounds z = θ[z_s2i[:]]
        @inbounds λ_all = θ[λ_all_s2i[:]]
        F_both!(result, z, x0, λ_all)
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
        @inbounds z = θ[z_s2i[:]]
        @inbounds λ_all = θ[λ_all_s2i[:]]
        J_vals!(data, z, x0, λ_all)
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

    @inbounds z = @view(θ[z_s2i[:]])

    (; status, info, θ, z, fres)
end