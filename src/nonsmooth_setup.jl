function gen_LP_data(A1::AbstractArray{T}, b1, A2, b2) where {T}
    m1 = length(b1)
    m2 = length(b2)
    A = [[A1; A2] ones(T, m1 + m2)]
    b = [b1; b2]
    q = [0, 0, 1.0]
    (A, b, q)
end

function g_col_single(xt, Ae, be, Ao, bo)
    sds = Dict()
    Aex, bex = shift_to(Ae, be, xt)
    AA, bb, qq = gen_LP_data(Aex, bex, Ao, bo)
    m1 = length(bex)
    m2 = length(bo)
    M = [zeros(Num, 3, 3) -AA'
        AA zeros(Num, m1 + m2, m1 + m2)]
    r = [qq; bb]
    all_active_inds = collect(1:m1+m2)
    Itr = powerset(all_active_inds) |> collect
    for active_inds in Itr
        length(active_inds) > 3 && continue
        assignment = [i ∈ active_inds for i in 1:m1+m2]
        try
            AA_active = collect(AA[active_inds, :])
            bb_active = collect(bb[active_inds])
            # TODO what if not unique primal? Need to resolve
            if length(active_inds) == 3
                zz = -AA_active \ bb_active
                #zz = -(AA_active'*AA_active)\AA_active'*bb_active
            else
                zz = -AA_active' * ((AA_active * AA_active') \ bb_active)
            end
            sd = zz[3]
            sds[active_inds] = sd
        catch err
            if err isa LinearAlgebra.SingularException
                continue
            else
                #@infiltrate
            end
        end
    end
    sds
end

function f_pack_single(xt, Ae, be, Q, q)
    fvals = Dict()
    Aex, bex = shift_to(Ae, be, xt)
    m = length(bex)
    M = [Q -Aex'
        Aex zeros(Num, m, m)]
    r = [q; bex]
    all_active_inds = collect(1:m)
    Itr = powerset(all_active_inds) |> collect
    for active_inds in Itr
        length(active_inds) > 2 && continue
        if length(active_inds) == 2 && abs(active_inds[1] - active_inds[2]) == 2
            continue
        end
        assignment = [i ∈ active_inds for i in 1:m]
        try
            AA_active = collect(Aex[active_inds, :])
            bb_active = collect(bex[active_inds])
            mm = length(active_inds)
            MM_active = [Q -AA_active'
                AA_active zeros(Num, mm, mm)]
            rr_active = [q; bb_active]
            rhs = -MM_active \ rr_active
            xx = rhs[1:2]
            fvals[active_inds] = 0.5 * xx' * Q * xx + xx' * q
        catch err
            if err isa LinearAlgebra.SingularException
                continue
            else
                #@infiltrate
            end
        end
    end
    fvals
end

function get_single_sd_ids(xt, Ae, be, Ao, bo, max_derivs)
    Aex, bex = shift_to(Ae, be, xt)
    AA, bb, qq = gen_LP_data(Aex, bex, Ao, bo)
    m1 = length(bex)
    m2 = length(bo)

    ret = solve_qp(UseOSQPSolver(); A=sparse(AA), l=-bb, q=qq, polish=true, verbose=false)
    #if ret.info.status_polish == -1
    #    @warn "not polished"
    #end
    primals = ret.x
    duals = -ret.y

    cons = AA * primals + bb
    I1 = duals .≥ 1e-2 .&& cons .< 1e-2
    I2 = duals .< 1e-2 .&& cons .< 1e-2
    I3 = duals .< 1e-2 .&& cons .≥ 1e-2

    all_inds = collect(1:m1+m2)
    weak_ind_options = powerset(all_inds[I2]) |> collect
    if length(weak_ind_options) > max_derivs
        @warn("More derivatives than accounted for! |I2| = $(sum(I2)). λ = $duals. cons = $cons")
        weak_ind_options = weak_ind_options[1:max_derivs]
    end
    if length(weak_ind_options) < max_derivs
        append!(weak_ind_options, [weak_ind_options[end] for _ in 1:max_derivs-length(weak_ind_options)])
    end

    assignments = map(weak_ind_options) do inds
        sort([all_inds[I1]; inds])
    end
end

function get_single_f_pack_ids(xt, Ae, be, Q, q, max_derivs, keys)
    Aex, bex = shift_to(Ae, be, xt)
    m = length(be)
    ret = solve_qp(UseOSQPSolver(); A=sparse(Aex), l=-bex, q, P=sparse(Q), polish=true, verbose=false)
    primals = ret.x
    duals = -ret.y
    cons = Aex * primals + bex
    I1 = duals .≥ 1e-2 .&& cons .< 1e-2
    I2 = duals .< 1e-2 .&& cons .< 1e-2
    I3 = duals .< 1e-2 .&& cons .≥ 1e-2

    all_inds = collect(1:m)

    if sum(I1) == 1
        for k in keys
            if all_inds[I1] ⊆ k
                for j in k
                    if !I1[j]
                        I2[j] = true
                    end
                end
            end
        end
    end

    weak_ind_options = powerset(all_inds[I2]) |> collect
    #if length(weak_ind_options) > max_derivs
    #    @warn("More derivatives than accounted for! |I2| = $(sum(I2)). λ = $duals. cons = $cons")
    #    weak_ind_options = weak_ind_options[1:max_derivs]
    #end

    assignments = map(weak_ind_options) do inds
        sort([all_inds[I1]; inds])
    end
    assignments = filter(a -> a ∈ keys, assignments)

    if length(assignments) > max_derivs
        #@infiltrate
        @warn("More derivatives than accounted for! |I2| = $(sum(I2)). λ = $duals. cons = $cons")
        assignments = assignments[1:max_derivs]
    end
    try
        if length(assignments) < max_derivs
            append!(assignments, [assignments[end] for _ in 1:max_derivs-length(assignments)])
        end
    catch e
        #@infiltrate
    end
    assignments
end


function setup_quick(ego_polys;
    T=50,
    dt=0.2,
    L=1.0,
    Q=0.01 * [1.0 0; 0 1],
    q=[0, 0.0],
    R=0.01 * I(3),
    p1_max=500.0,
    p2_min=-500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=π / 4,
    sides_per_poly=4,
    derivs_per_sd=4,
    derivs_per_fv=4,
    N_polys=4,
    rng=MersenneTwister(420))

    #Ae = Symbolics.@variables(Ae[1:sides_per_poly,1:2])[1] |> Symbolics.scalarize
    #be = Symbolics.@variables(be[1:sides_per_poly])[1] |> Symbolics.scalarize

    N_ego_polys = length(ego_polys)
    Ao = Symbolics.@variables(Ao[1:sides_per_poly, 1:2])[1] |> Symbolics.scalarize
    bo = Symbolics.@variables(bo[1:sides_per_poly])[1] |> Symbolics.scalarize

    z = Symbolics.@variables(z[1:9*T])[1] |> Symbolics.scalarize
    xt = Symbolics.@variables(xt[1:6])[1] |> Symbolics.scalarize
    λsd = Symbolics.@variables(λsd)
    α = Symbolics.@variables(α)
    β = Symbolics.@variables(β)
    x0 = Symbolics.@variables(x0[1:6])[1] |> Symbolics.scalarize

    sds = map(ego_polys) do P
        Ae = collect(P.A)
        be = P.b
        g_col_single(xt, Ae, be, Ao, bo)
    end

    num_sd_cons = T * N_polys * N_ego_polys
    num_sd_mults = T * N_polys * derivs_per_sd * N_ego_polys

    fvals = map(ego_polys) do P
        Ae = collect(P.A)
        be = P.b
        f_pack_single(xt, Ae, be, Q, q)
    end
    fkeys = map(fvals) do fv
        collect(keys(fv))
    end
    num_f_mults = derivs_per_fv * N_ego_polys

    cost_nom = f(z, T, R)
    cons_dyn = g_dyn(z, x0, T, dt, L)
    cons_env = g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)
    cons_nom = [cons_dyn; cons_env]

    α_f = Symbolics.@variables(α_f[1:num_f_mults])[1] |> Symbolics.scalarize
    β_sd = Symbolics.@variables(β_sd[1:num_sd_mults])[1] |> Symbolics.scalarize
    slacks = Symbolics.@variables(slacks[1:num_sd_cons+N_ego_polys])[1] |> Symbolics.scalarize
    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize
    λ_col = Symbolics.@variables(λ_col[1:num_sd_cons])[1] |> Symbolics.scalarize

    θ = [z; α_f; β_sd; slacks; λ_nom; λ_col]

    lag = cost_nom - cons_nom' * λ_nom
    grad_lag = Symbolics.gradient(lag, z)

    simplex_cons = Num[]
    for i in 1:N_ego_polys
        push!(simplex_cons, 1.0 - sum(α_f[(i-1)*derivs_per_fv+1:i*derivs_per_fv]))
    end
    for i in 1:N_ego_polys
        offset = (i - 1) * T * N_polys * derivs_per_sd
        for t in 1:T
            for j in 1:N_polys
                push!(simplex_cons, 1.0 - sum(β_sd[(1:derivs_per_sd).+(offset+(t-1)*N_polys*derivs_per_sd+(j-1)*derivs_per_sd)]))
            end
        end
    end
    F_nom = [grad_lag; zeros(Num, num_f_mults + num_sd_mults); simplex_cons; cons_nom; zeros(Num, num_sd_cons)]

    l = [fill(-Inf, length(grad_lag)); fill(0.0, num_f_mults + num_sd_mults); fill(-Inf, length(simplex_cons)); fill(-Inf, length(cons_dyn)); zeros(length(cons_env)); zeros(num_sd_cons)]
    u = [fill(+Inf, length(grad_lag)); fill(Inf, num_f_mults + num_sd_mults); fill(Inf, length(simplex_cons)); fill(Inf, length(cons_dyn)); fill(Inf, length(cons_env)); fill(Inf, num_sd_cons)]
    n = length(l)

    Jnom = Symbolics.sparsejacobian(F_nom, θ)
    (Jnom_rows, Jnom_cols, Jnom_vals) = findnz(Jnom)
    Jnom_buf = zeros(length(Jnom_vals))
    get_Fnom! = Symbolics.build_function(F_nom, z, x0, λ_nom, α_f, β_sd; expression=Val(false))[2]
    get_Jnom_vals = Symbolics.build_function(Jnom_vals, z, x0, λ_nom, α_f, β_sd; expression=Val(false))[2]

    sd_grad_xt_lag_funcs = []
    sd_funcs = []

    get_lag = Dict()
    get_sd = Dict()
    get_Jlag = Dict()
    get_Jsd = Dict()

    get_gfv = Dict()
    get_Jfv = Dict()

    println(length(sds))
    @info "Generating symbolic solutions to sd calculation"
    for (i, sds_i) in enumerate(sds)
        @showprogress for (k, sd) in sds_i
            lag = Symbolics.gradient(-sd * λsd[1] * β[1], xt; simplify=false)
            get_lag[i, k] = Symbolics.build_function(lag, xt, Ao, bo, β, λsd; expression=Val(false))[2]
            get_sd[i, k] = Symbolics.build_function(β[1] * sd, xt, Ao, bo, β, λsd; expression=Val(false))

            Jlag = Symbolics.sparsejacobian(lag, [xt; β; λsd]; simplify=false)
            Jsd = Symbolics.sparsejacobian([β[1] * sd], [xt; β]; simplify=false)

            Jlag_rows, Jlag_cols, Jlag_vals = findnz(Jlag)
            if length(Jlag_vals) == 0
                Jlag_rows = [1]
                Jlag_cols = [1]
                Jlag_vals = Num[0.0]
            end
            Jsd_rows, Jsd_cols, Jsd_vals = findnz(Jsd)
            if length(Jsd_vals) == 0
                Jsd_rows = [1]
                Jsd_cols = [1]
                Jsd_vals = Num[0.0]
            end
            get_Jlag[i, k] = (Jlag_rows, Jlag_cols, Symbolics.build_function(Jlag_vals, xt, Ao, bo, β, λsd; expression=Val(false))[2], zeros(length(Jlag_rows)))
            get_Jsd[i, k] = (Jsd_rows, Jsd_cols, Symbolics.build_function(Jsd_vals, xt, Ao, bo, β, λsd; expression=Val(false))[2], zeros(length(Jsd_rows)))
        end
    end
    @info "Generating symbolic solutions to f calculation"
    for (i, fvals_i) in enumerate(fvals)
        @showprogress for (k, fv) in fvals_i
            gfv = Symbolics.gradient(fv * α[1], xt; simplify=false)
            get_gfv[i, k] = Symbolics.build_function(gfv, xt, α; expression=Val(false))[2]

            Jfv = Symbolics.sparsejacobian(gfv, [xt; α]; simplify=false)
            Jfv_rows, Jfv_cols, Jfv_vals = findnz(Jfv)
            get_Jfv[i, k] = (Jfv_rows, Jfv_cols, Symbolics.build_function(Jfv_vals, xt, α; expression=Val(false))[2], zeros(length(Jfv_rows)))
        end
    end

    Aes = [deepcopy(P.A) for P in ego_polys]
    bes = [deepcopy(P.b) for P in ego_polys]

    col_offset = length(z) + length(α_f) + length(β_sd) + length(slacks) + length(λ_nom)
    β_offset = length(z) + length(α_f)
    lag_buf = zeros(6)
    function fill_F!(F, z, x0, polys, α_f, β_sd, λ_nom, λ_col)
        F .= 0.0
        get_Fnom!(F, z, x0, λ_nom, α_f, β_sd)
        for i in 1:N_ego_polys
            for t in 1:T
                xt_inds = (t-1)*9+1:(t-1)*9+6
                @inbounds xt = z[xt_inds]
                for (e, P) in enumerate(polys)
                    Ao = P.A
                    bo = P.b
                    λ_ind = (i - 1) * T * N_polys + (t - 1) * N_polys + e
                    β_inds = (1:derivs_per_sd) .+ (((i - 1) * T * N_polys * derivs_per_sd) + (t - 1) * N_polys * derivs_per_sd + (e - 1) * derivs_per_sd)
                    @inbounds λte = λ_col[λ_ind]
                    @inbounds βte = β_sd[β_inds]
                    assignments = get_single_sd_ids(xt, Aes[i], bes[i], Ao, bo, derivs_per_sd)
                    for (ee, assignment) in enumerate(assignments)
                        if length(assignment) > 3
                            assignment = assignment[1:3]
                        end
                        get_lag[i, assignment](lag_buf, xt, Ao, bo, βte[ee], λte)
                        #@infiltrate any(isnan.(βte[ee]))
                        βsd = get_sd[i, assignment](xt, Ao, bo, βte[ee], λte)

                        #if any(isnan.(lag_buf))
                        #    # lag_buf can be nan
                        #    # needs permanent solution
                        #    #@infiltrate
                        #    xt[3] += 1e-15;
                        #    get_lag[i, assignment](lag_buf, xt, Ao, bo, βte[ee], λte)
                        #    @infiltrate any(isnan.(lag_buf))
                        #end

                        F[xt_inds] .+= lag_buf
                        F[λ_ind+col_offset] += βsd
                        #F[λte_inds .+ offset] .+= sd*λ_normalized[ee]
                    end
                end

                if t == T
                    assignments_f = get_single_f_pack_ids(xt, Aes[i], bes[i], Q, q, derivs_per_fv, fkeys[i])
                    α_inds = (i-1)*derivs_per_fv+1:i*derivs_per_fv
                    @inbounds αi = α_f[α_inds]
                    for (ee, assignment) in enumerate(assignments_f)
                        get_gfv[i, assignment](lag_buf, xt, αi[ee])
                        F[xt_inds] .+= lag_buf
                    end
                end
            end
        end
        nothing
    end

    function get_J_both(JJ, z, x0, polys, α_f, β_sd, λ_nom, λ_col)
        JJ.nzval .= 1e-16
        get_Jnom_vals(Jnom_buf, z, x0, λ_nom, α_f, β_sd)
        JJ .+= sparse(Jnom_rows, Jnom_cols, Jnom_buf, n, n)

        for i in 1:N_ego_polys
            for t in 1:T
                xt_inds = (t-1)*9+1:(t-1)*9+6
                @inbounds xt = z[xt_inds]
                for (e, P) in enumerate(polys)
                    Ao = collect(P.A)
                    bo = P.b
                    λ_ind = (i - 1) * T * N_polys + (t - 1) * N_polys + e
                    β_inds = (1:derivs_per_sd) .+ (((i - 1) * T * N_polys * derivs_per_sd) + (t - 1) * N_polys * derivs_per_sd + (e - 1) * derivs_per_sd)
                    @inbounds λte = λ_col[λ_ind]
                    @inbounds βte = β_sd[β_inds]
                    #λte_inds = (1:derivs_per_sd) .+ ((t-1)*N_polys*derivs_per_sd+(e-1)*derivs_per_sd)
                    #@inbounds λte = λ_col[λte_inds]
                    assignments = get_single_sd_ids(xt, Aes[i], bes[i], Ao, bo, derivs_per_sd)
                    for (ee, assignment) in enumerate(assignments)
                        if length(assignment) > 3
                            assignment = assignment[1:3]
                        end
                        #Jlag = get_Jlag[assignment](xt,Ae,be,Ao,bo,λte[ee])
                        #Jsd = get_Jsd[assignment](xt,Ae,be,Ao,bo,λte[ee])
                        Jlag_rows, Jlag_cols, Jlag_vals, Jlag_buf = get_Jlag[i, assignment]
                        Jsd_rows, Jsd_cols, Jsd_vals, Jsd_buf = get_Jsd[i, assignment]
                        Jlag_buf .= 0.0
                        Jsd_buf .= 0.0
                        Jlag_vals(Jlag_buf, xt, Ao, bo, βte[ee], λte)
                        Jsd_vals(Jsd_buf, xt, Ao, bo, βte[ee], λte)
                        Jlag = sparse(Jlag_rows, Jlag_cols, Jlag_buf, 6, 8)
                        Jsd = sparse(Jsd_rows, Jsd_cols, Jsd_buf, 1, 7)

                        #@infiltrate any(isnan.(Jlag))

                        #if any(isnan.(Jlag_buf))
                        #    # Jlag and Jsd can be nan
                        #    # needs permanent solution
                        #    @infiltrate
                        #    xt[3] += 1e-15;
                        #    Jlag_vals(Jlag_buf, xt, Ao, bo, βte[ee], λte)
                        #    Jsd_vals(Jsd_buf, xt, Ao, bo, βte[ee], λte)
                        #    @infiltrate any(isnan.(Jlag_buf))
                        #    @infiltrate any(isnan.(Jsd_buf))
                        #end

                        @inbounds JJ[xt_inds, xt_inds] .+= Jlag[1:6, 1:6]
                        @inbounds JJ[xt_inds, β_inds[ee]+β_offset] .+= Jlag[1:6, 7]
                        @inbounds JJ[xt_inds, λ_ind+col_offset] .+= Jlag[1:6, 8]
                        @inbounds JJ[λ_ind+col_offset, xt_inds] .+= Jsd[1:6]
                        @inbounds JJ[λ_ind+col_offset, β_inds[ee]+β_offset] += Jsd[7]
                    end
                end
                if t == T
                    assignments_f = get_single_f_pack_ids(xt, Aes[i], bes[i], Q, q, derivs_per_fv, fkeys[i])
                    α_inds = (i-1)*derivs_per_fv+1:i*derivs_per_fv
                    @inbounds αi = α_f[α_inds]
                    for (ee, assignment) in enumerate(assignments_f)
                        Jfv_rows, Jfv_cols, Jfv_vals, Jfv_buf = get_Jfv[i, assignment]
                        Jfv_vals(Jfv_buf, xt, αi[ee])
                        Jfv = sparse(Jfv_rows, Jfv_cols, Jfv_buf, 6, 7)
                        #@infiltrate any(isnan.(Jfv))

                        @inbounds JJ[xt_inds, xt_inds] .+= Jfv[:, 1:length(xt_inds)]
                        @inbounds JJ[xt_inds, length(z)+α_inds[ee]] .+= Jfv[:, end]
                    end
                end
            end
        end
        nothing
    end

    J_example = sparse(Jnom_rows, Jnom_cols, ones(length(Jnom_cols)), n, n)
    for i in 1:N_ego_polys
        for t in 1:T
            xt_inds = (t-1)*9+1:(t-1)*9+6
            for e in 1:N_polys
                λ_ind = (i - 1) * T * N_polys + (t - 1) * N_polys + e
                β_inds = (1:derivs_per_sd) .+ (((i - 1) * T * N_polys * derivs_per_sd) + (t - 1) * N_polys * derivs_per_sd + (e - 1) * derivs_per_sd)
                J_example[xt_inds, xt_inds] .= 1.0
                J_example[xt_inds, col_offset+λ_ind] .= 1.0
                J_example[col_offset+λ_ind, xt_inds] .= 1.0
                for l in β_inds
                    J_example[xt_inds, β_offset+l] .= 1.0
                    J_example[col_offset+λ_ind, β_offset+l] = 1.0
                end
            end
            if t == T
                J_example[xt_inds, length(z)+1:length(z)+length(α_f)] .= 1.0
            end
        end
    end

    @info "Forcing compilation"
    n_z = length(z)
    n_nom = length(λ_nom)
    n_α = length(α_f)
    n_β = length(β_sd)
    n_s = length(slacks)
    n_col = length(λ_col)
    xtr = randn(6)
    Aer = randn(sides_per_poly, 2)
    ber = randn(sides_per_poly)
    Aor = randn(sides_per_poly, 2)
    bor = randn(sides_per_poly)
    lag_buf = zeros(6)
    λsdr = randn()
    αr = randn()
    βr = randn()
    for i in 1:N_ego_polys
        @showprogress for k in collect(keys(sds[i]))
            get_lag[i, k](lag_buf, xtr, Aor, bor, βr, λsdr)
            ss = get_sd[i, k](xtr, Aor, bor, βr, λsdr)
            get_Jlag[i, k][3](get_Jlag[i, k][4], xtr, Aor, bor, βr, λsdr)
            get_Jsd[i, k][3](get_Jsd[i, k][4], xtr, Aor, bor, βr, λsdr)
            #ll = get_lag[k](xtr,Aer,ber,Aor,bor, λsdr)
            #ss = get_sd[k](xtr,Aer,ber,Aor,bor, λsdr)
            #J1 = get_Jlag[k](xtr,Aer,ber,Aor,bor,λsdr)
            #J2 = get_Jsd[k](xtr,Aer,ber,Aor,bor,λsdr)
        end
        @showprogress for k in collect(keys(fvals[i]))
            gfv = get_gfv[i, k](lag_buf, xtr, αr)
            Jfv = get_Jfv[i, k][3](get_Jfv[i, k][4], xtr, αr)
        end
    end

    #@infiltrate

    return (; fill_F!,
        get_J_both,
        J_example,
        l,
        u,
        T,
        n_z,
        n_α,
        n_β,
        n_s,
        n_nom,
        n_col,
        N_polys,
        ego_polys,
        sides_per_poly,
        p1_max,
        p2_min)
end



function solve_quick(prob, x0, polys; θ0=nothing)
    (; fill_F!, get_J_both, J_example, ego_polys, l, u, T, n_z, n_α, n_β, n_s, n_nom, n_col, N_polys, sides_per_poly, p1_max, p2_min) = prob

    @assert length(polys) == N_polys

    n = length(l)
    @assert n == n_z + n_α + n_β + n_s + n_nom + n_col

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())


    xxts = Dict()
    #for t in 10:10:T
    for i in 1:length(ego_polys)
        xx = x0[1:3]
        Aeb = shift_to(ego_polys[i].A, ego_polys[i].b, xx)
        self_poly = ConvexPolygon2D(Aeb[1], Aeb[2])
        plot!(ax, self_poly; color=:blue)
        for t in 5:1:T-1
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

    #colors = [:red, :orange, :yellow, :green]
    colors = [:red for _ in 1:N_polys]
    for (P, c) in zip(polys, colors)
        plot!(ax, P; color=c)
    end

    #lines!(ax, [-p1_max, p1_max],[p2_min,p2_min], color=:black)
    #lines!(ax, [-p1_max, -p1_max], [p2_min, 10], color=:black)
    #lines!(ax, [p1_max, p1_max], [p2_min, 10], color=:black)
    display(fig)

    if isnothing(θ0)
        θ0 = zeros(n)
        derivs_per_fv = Int(n_α / length(ego_polys))
        derivs_per_sd = Int(n_β / (T * length(ego_polys) * N_polys))
        for t in 1:T
            θ0[(t-1)*9+1:(t-1)*9+6] = x0
        end
        θ0[n_z+1:n_z+n_α] .= 1.0 / derivs_per_fv
        θ0[n_z+n_α+1:n_z+n_α+n_β] .= 1.0 / derivs_per_sd
    end

    JJ = J_example
    J_col = JJ.colptr[1:end-1]
    J_len = diff(JJ.colptr)
    J_row = JJ.rowval
    nnz_total = length(JJ.nzval)

    function F(n, θ, result)
        result .= 0.0
        #@inbounds z = @view(θ[1:n_z])
        #@inbounds λ_nom = @view(θ[n_z+1:n_z+n_nom])
        #@inbounds λ_col = @view(θ[n_z+n_nom+1:n_z+n_nom+n_col])
        @inbounds z = θ[1:n_z]
        @inbounds α_f = θ[n_z+1:n_z+n_α]
        @inbounds β_sd = θ[n_z+n_α+1:n_z+n_α+n_β]
        @inbounds λ_nom = θ[n_z+n_α+n_β+n_s+1:n_z+n_α+n_β+n_s+n_nom]
        @inbounds λ_col = θ[n_z+n_α+n_β+n_s+n_nom+1:n_z+n_α+n_β+n_s+n_nom+n_col]
        fill_F!(result, z, x0, polys, α_f, β_sd, λ_nom, λ_col)
        #for t in 10:10:T
        for i in 1:length(ego_polys)
            for t in 5:5:T
                xxts[i, t][] = copy(θ[(t-1)*9+1:(t-1)*9+6])
            end
        end
        Cint(0)
    end
    function J(n, nnz, θ, col, len, row, data)
        @assert nnz == nnz_total
        data .= 0.0
        @inbounds z = θ[1:n_z]
        @inbounds α_f = θ[n_z+1:n_z+n_α]
        @inbounds β_sd = θ[n_z+n_α+1:n_z+n_α+n_β]
        @inbounds λ_nom = θ[n_z+n_α+n_β+n_s+1:n_z+n_α+n_β+n_s+n_nom]
        @inbounds λ_col = θ[n_z+n_α+n_β+n_s+n_nom+1:n_z+n_α+n_β+n_s+n_nom+n_col]
        #@infiltrate any(isnan.(z))
        get_J_both(JJ, z, x0, polys, α_f, β_sd, λ_nom, λ_col)
        #@infiltrate any(isnan.(JJ))
        col .= J_col
        len .= J_len
        row .= J_row
        data .= JJ.nzval
        Cint(0)
    end


    buf = zeros(n)
    buf2 = zeros(n)
    Jbuf = zeros(nnz_total)

    w = randn(length(θ0))
    w = copy(θ0)

    F(n, w, buf)
    J(n, nnz_total, w, zero(J_col), zero(J_len), zero(J_row), Jbuf)

    Jrows, Jcols, _ = findnz(J_example)

    Jnum = sparse(Jrows, Jcols, Jbuf)
    Jnum2 = spzeros(n, n)
    #@info "Testing Jacobian accuracy numerically"
    #@showprogress for ni in 1:n
    #    wi = copy(w)
    #    wi[ni] += 1e-5
    #    F(n, wi, buf2)
    #    Jnum2[:,ni] = sparse((buf2-buf) ./ 1e-5)
    #end
    #@info "Jacobian error is $(norm(Jnum2-Jnum))"

    #@infiltrate
    #Main.@infiltrate

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
        restart_limit=0,
        jacobian_data_contiguous=true,
        cumulative_iteration_limit=50_000,
        convergence_tolerance=5e-4
        #proximal_perturbation=0,
        #crash_method="pnewton",
        #crash_nbchange_limit=10,
        #nms_initial_reference_factor=2,
        #crash_searchtype="arc",
        #nms_searchtype="arc",
        #gradient_searchtype="arc",
        #lemke_search_type="slack"
    )

    fres = zeros(n)

    F(n, θ, fres)

    display(fig)

    @inbounds z = @view(θ[1:n_z])

    (; status, info, θ, z, fres)
end




