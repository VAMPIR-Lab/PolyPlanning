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
    T=1,
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
    z_s2i, dyn_cons_s2i, env_cons_s2i, sd_cons_s2i = get_sub2idxs((n_xu, T), (n_dyn_cons), (n_env_cons), (n_sd_cons_per_col, n_ego, n_obs, T))

    z_sym = Symbolics.@variables(z[1:n_z])[1] |> Symbolics.scalarize
    x0_sym = Symbolics.@variables(x0[1:n_x])[1] |> Symbolics.scalarize
    xt_sym = Symbolics.@variables(xt[1:n_x])[1] |> Symbolics.scalarize

    cost = f(z_sym, T, R_cost, Q_cost)
    dyn_cons = g_dyn(z_sym, x0_sym, T, dt)
    env_cons = g_env(z_sym, T, p1_max, p2_min, u1_max, u2_max, u3_max)

    # indexing consistency
    @assert length(z_s2i) == n_z
    @assert length(dyn_cons_s2i) == length(dyn_cons)
    @assert length(env_cons_s2i) == length(env_cons)
    @assert length(sd_cons_s2i) == n_sd_cons

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
            g_col_single(xt_sym, A_ego, b_ego, centr_ego, A_obs, b_obs, centr_obs)
        end
    end

    @infiltrate
    all_cons = [dyn_cons; env_cons;]


    # dynamic and environment constraints
    λ_nom = Symbolics.@variables(λ_nom[1:length(all_cons)])[1] |> Symbolics.scalarize
    # collision avoidance constraints
    λ_col = Symbolics.@variables(λ_col[1:num_sd_cons])[1] |> Symbolics.scalarize

    θ = [z; α_f; β_sd; slacks; λ_nom; λ_col]

    for i in 1:n_ego
        offset = (i - 1) * T * n_obs * derivs_per_sd
        for t in 1:T
            for j in 1:n_obs
                push!(simplex_cons, 1.0 - sum(β_sd[(1:derivs_per_sd).+(offset+(t-1)*n_obs*derivs_per_sd+(j-1)*derivs_per_sd)]))
            end
        end
    end
    F_nom = [grad_lag; zeros(Num, num_f_mults + num_sd_mults); simplex_cons; all_cons; zeros(Num, num_sd_cons)]

    l = [fill(-Inf, length(grad_lag)); fill(0.0, num_f_mults + num_sd_mults); fill(-Inf, length(simplex_cons)); fill(-Inf, length(dyn_cons)); zeros(length(env_cons)); zeros(num_sd_cons)]
    u = [fill(+Inf, length(grad_lag)); fill(Inf, num_f_mults + num_sd_mults); fill(Inf, length(simplex_cons)); fill(Inf, length(dyn_cons)); fill(Inf, length(env_cons)); fill(Inf, num_sd_cons)]
    n = length(l)

    Jnom = Symbolics.sparsejacobian(F_nom, θ)

    (Jnom_rows, Jnom_cols, Jnom_vals) = findnz(Jnom)
    Jnom_buf = zeros(length(Jnom_vals))
    get_Fnom! = Symbolics.build_function(F_nom, z, x0, λ_nom, α_f, β_sd; expression=Val(false))[2]
    get_Jnom_vals = Symbolics.build_function(Jnom_vals, z, x0, λ_nom, α_f, β_sd; expression=Val(false))[2]

    get_lag = Dict()
    get_sd = Dict()
    get_Jlag = Dict()
    get_Jsd = Dict()

end

