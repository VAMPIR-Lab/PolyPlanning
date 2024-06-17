#function gen_LP_data(xt, Ae::AbstractArray{T}, be, centroide, Ao, bo, centroido; is_newsd=false) where {T}
#    me = length(be)
#    mo = length(bo)
#    A = [[Ae; Ao] ones(T, me + mo)]
#    b = [be; bo]
#    q = [0, 0, 1.0]
#    # print("old")
#    (A, b, q)
#end

function gen_LP_data(xt, Ae::AbstractArray{T}, be, centroide, Ao, bo, centroido; is_newsd=false) where {T}
    A = [Ae Ae*centroide+be
        Ao Ao*centroido+bo]
    b = [be; bo]
    q = [0, 0, 1.0]
    # print("new")
    (A, b, q)
end

# function gen_LP_data(xt, Ae::AbstractArray{T}, be, centroide, Ao, bo, centroido; is_newsd=false) where {T}
#     if !is_newsd
#         me = length(be)
#         mo = length(bo)
#         A = [[Ae; Ao] ones(T, me + mo)]
#         b = [be; bo]
#         q = [0, 0, 1.0]
#     else
#         A = [Ae Ae*centroide+be;
#             Ao Ao*centroido+bo]
#         b = [be; bo]
#         q = [0, 0, 1.0]
#     end
#     (A, b, q)
# end

function g_col_single(xt, Ae, be, centroide, Ao, bo, centroido; is_newsd=false)
    sds = Dict()
    Aex, bex = shift_to(Ae, be, xt)
    R = [cos(xt[3]) sin(xt[3])
        -sin(xt[3]) cos(xt[3])]
    centroidex = xt[1:2] + R * centroide
    AA, bb, qq = gen_LP_data(xt, Aex, bex, centroidex, Ao, bo, centroido; is_newsd=is_newsd)
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
                #if linear system is underdetermined, use minimum norm solution (calculated by right inverse)
                #Note: every solution has the same sd, i.e., zz[3], but different zz[1:2]
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
                @warn(err)
            end
        end
    end
    fvals
end

function plot_inflated(xt, Ae, be, centroide, Ao, bo, centroido, p_intercept, sd; fig=Figure(), ax=Axis(fig[1, 1], aspect=DataAspect()), θ=[], is_displaying=true)

    # draw polys
    Aeb = shift_to(Ae, be, xt)
    Aex = Aeb[1]
    bex = Aeb[2]
    ego_shifted = ConvexPolygon2D(Aex, bex)
    plot!(ax, ego_shifted; color=:blue)
    obs = ConvexPolygon2D(Ao, bo)
    plot!(ax, obs; color=:red)

    # draw inflated
    be_inflated = bex + sd * (Aex * centroide + bex)
    bo_inflated = bo + sd * (Ao * centroido + bo)
    ego_inflated = ConvexPolygon2D(Aex, be_inflated)
    obs_inflated = ConvexPolygon2D(Ao, bo_inflated)

    plot!(ax, ego_inflated; color=:lightblue, linewidth=2)
    plot!(ax, obs_inflated; color=:pink, linewidth=2)

    # draw the intercept and centroids
    scatter!(ax, centroide[1], centroide[2]; color=:blue)
    scatter!(ax, centroido[1], centroido[2]; color=:red)
    scatter!(ax, p_intercept[1], p_intercept[2]; color=:green)

    if is_displaying
        display(fig)
    end

    (fig, ax)
end

# some index combinations is not feasible, i.e. all indices are from m1 or m2
function get_single_sd_ids(xt, Ae, be, centroide, Ao, bo, centroido, max_derivs, t; is_newsd=false)
    Aex, bex = shift_to(Ae, be, xt)
    R = [cos(xt[3]) sin(xt[3])
        -sin(xt[3]) cos(xt[3])]
    centroidex = xt[1:2] + R * centroide
    AA, bb, qq = gen_LP_data(xt, Aex, bex, centroidex, Ao, bo, centroido; is_newsd=is_newsd)
    m1 = length(bex)
    m2 = length(bo)

    # use OSQP solver
    #ret = solve_qp(UseOSQPSolver(); A=sparse(AA), l=-bb, q=qq, polish=true, verbose=false, max_iter=4e4, eps_abs=1e-5)#, eps_abs=1e-5, eps_rel=1e-5, max_iter=1e6)

    #if ret.info.status_val != 1
    #    @warn ret.info.status
    #    Main.@infiltrate
    #    #plot_inflated(xt, Ae, be, centroidex, Ao, bo, centroido, ret.x[1:2], ret.x[3])
    #end
    #if ret.info.status_polish == -1
    #    @warn "not polished"
    #end
    #primals = ret.x
    #duals = -ret.y

    #cons = AA * primals + bb
    #I1 = duals .≥ 1e-2 .&& cons .< 1e-2
    #I2 = duals .< 1e-2 .&& cons .< 1e-2
    #I3 = duals .< 1e-2 .&& cons .≥ 1e-2

    ## use JuMP and Clp solver
    model = JuMP.Model(Clp.Optimizer)
    JuMP.set_attribute(model, "LogLevel", 0) # disable printing log
    JuMP.@variable(model, xx[1:3])
    JuMP.@constraint(model, constraint, AA * xx + bb .>= 0)
    JuMP.@objective(model, Min, qq' * xx)
    JuMP.optimize!(model)

    status = JuMP.termination_status(model)
    if status != OPTIMAL
        duals = zeros(m1 + m2)
        cons = duals
        xx = JuMP.value.(xx)
        @warn status
        #Main.@infiltrate
        #plot_inflated(xt, Ae, be, centroidex, Ao, bo, centroido, xx[1:2], xx[3])
    else
        primals = JuMP.value.(xx)
        duals = JuMP.dual.(constraint)
        cons = AA * primals + bb
    end

    # old
    ####
    # primal and dual tolerance is 1e-7
    #I1 = duals .≥ 1e-6 .&& cons .< 1e-6
    #I2 = duals .< 1e-6 .&& cons .< 1e-6
    #I3 = duals .< 1e-6 .&& cons .≥ 1e-6

    #all_inds = collect(1:m1+m2)
    #weak_ind_options = powerset(all_inds[I2]) |> collect
    ## the first set is an empty set, and is always preserved
    #if length(weak_ind_options) > max_derivs
    #    @warn("More sd derivatives than accounted for! |I2| = $(sum(I2)). λ = $duals. cons = $cons")
    #    weak_ind_options = weak_ind_options[1:max_derivs]
    #end
    ## what if weak_ind_options is empty?
    #if length(weak_ind_options) < max_derivs
    #    append!(weak_ind_options, [weak_ind_options[end] for _ in 1:max_derivs-length(weak_ind_options)])
    #end

    ##indices of I1 + subset of I2
    #assignments = map(weak_ind_options) do inds
    #    sort([all_inds[I1]; inds])
    #end

    # new 
    #####
    assignments = Set{Set{Int}}()
    verts = []
    explored_verts = []
    candidates = []

    push!(candidates, [primals; duals])

    #@infiltrate

    while !isempty(candidates)
        #@infiltrate
        v = pop!(candidates)
        pp = v[1:3]
        dd = v[4:end]
        push!(explored_verts, v)
        cons = AA * pp + bb
        ### primal and dual tolerance is 1e-7
        tol = 1e-3
        I1 = dd .≥ tol .&& cons .< tol
        I2 = dd .< tol .&& cons .< tol
        I3 = dd .< tol .&& cons .≥ tol

        all_inds = collect(1:m1+m2)
        weak_ind_options = powerset(all_inds[I2]) |> collect
        # the first set is an empty set, and is always preserved
        if length(weak_ind_options) > max_derivs
            #@warn("More sd derivatives than accounted for! |I2| = $(sum(I2)). λ = $duals. cons = $cons")
            weak_ind_options = weak_ind_options[1:max_derivs]
        end
        # what if weak_ind_options is empty?
        if length(weak_ind_options) < max_derivs
            append!(weak_ind_options, [weak_ind_options[end] for _ in 1:max_derivs-length(weak_ind_options)])
        end

        for inds in weak_ind_options
            ass = Set{Int}([all_inds[I1]; inds])
            ma = length(ass) # number of active indices
            ass_v = ass |> collect |> sort
            # B_eq * [primals; duals_active] + b_eq = 0
            # B_ineq * duals_inactive + b_ineq ≥ 0
            B_eq = [spzeros(3, 3) -AA[ass_v, :]'
                AA[ass_v, :] spzeros(ma, ma)]
            not_ass = setdiff(all_inds, ass) # is this always ordered
            B_ineq = [AA[not_ass, :] spzeros(m1 + m2 - ma, ma)]
            b_eq = [qq; bb[ass_v]]
            b_ineq = bb[not_ass]

            #if !(B_eq * [pp; dd[ass_v]] + b_eq |> norm < 1e-5)
            #    @infiltrate
            #end
            #@assert B_eq * [pp; dd[ass_v]] + b_eq |> norm < 1e-5
            #@assert all(B_ineq * [pp; dd[ass_v]] + b_ineq .≥ -1e-5)

            #@infiltrate
            # B_eq2 * [primals; duals_active; duals_inactive] + b_eq = 0
            # B_ineq2 * [primals; duals_active; duals_inactive] + b_ineq ≥ 0 
            B_eq2 = [B_eq zeros(3 + ma, m1 + m2 - ma)]
            B_ineq2 = [B_ineq zeros(m1 + m2 - ma, m1 + m2 - ma)]
            B_big = [B_eq2; B_ineq2]

            # hyperplanes, i.e., equality constraints
            hps = map(1:ma+3) do i
                HyperPlane(-B_eq2[i, :], b_eq[i])
            end
            
            # halfspaces, i.e., inequality constraints
            hss = map(1:m1+m2-ma) do i
                HalfSpace(-B_ineq2[i, :], b_ineq[i])
            end
            hr = hrep(hps, hss) # H-representation for polyhedron 
            poly = Polyhedra.polyhedron(hr)
            vr = vrep(poly) # V-representation

            local_verts = Polyhedra.points(vr)
            @infiltrate isempty(local_verts)

            # need to consider the order of indices
            for v in local_verts
                #@infiltrate
                xv = v[1:3]
                λv_eq = v[4:end]
                λv = zeros(m1 + m2)
                λv[ass_v] = λv_eq[1:ma]
                vv = [xv; λv]
                if !any(isapprox(vv, vi; atol=1e-6) for vi in [explored_verts; verts])
                    push!(verts, vv)
                    push!(candidates, vv)
                end
            end

            if length(ass) ≥ 3
                push!(assignments, ass)
            end
        end
    end

    #@infiltrate





    #function does_it_exist(as, v)
    #    for (i, a) in enumerate(as)
    #        if a == v
    #            return i
    #        end
    #    end
    #    return 0
    #end

    #new_assignments = copy(assignments)

    ##for (i, a) in enumerate(assignments)
    #i = does_it_exist(assignments, [1, 8])
    #if i > 0
    #    # if we have [1,2, 8] already...
    #    if does_it_exist(assignments, [1, 2, 8]) == 0
    #        new_assignments[i] = [1, 2, 8]
    #        # if we have [1,4,8] already
    #    elseif does_it_exist(assignments, [1, 4, 8]) == 0
    #        new_assignments[i] = [1, 4, 8]
    #    end
    #end
    #end


    #if t == 20
    #plot_inflated(xt, Ae, be, centroidex, Ao, bo, centroido, primals[1:2], primals[3])
    #@info "old:"
    #@info assignments
    #@info "hacked:"
    #@info new_assignments
    #Main.@infiltrate
    #end


    v_assignments = []
    for a in assignments
        push!(v_assignments, sort(collect(a)))
    end
    #@infiltrate

    if length(v_assignments) < 1
        @infiltrate
    end
    v_assignments


    # do we need this??
    #if length(v_assignments) > max_derivs
    #    v_assignments[1:max_derivs]
    #end
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
        @warn("More f_pack derivatives than accounted for! |I2| = $(sum(I2)). λ = $duals. cons = $cons")
        assignments = assignments[1:max_derivs]
    end
    # if length(assignments)=0, nothing is added
    try
        if length(assignments) < max_derivs
            append!(assignments, [assignments[end] for _ in 1:max_derivs-length(assignments)])
        end
    catch e
        warn(e)
    end
    assignments
end


function setup_quick(ego_polys;
    T=50,
    dt=0.2,
    Q=0.01 * [1.0 0; 0 1],
    q=[0, 0.0],
    enable_fvals=false,
    Rf=1e-3 * I(3),
    Qf=1e-3 * I(2),
    p1_max=500.0,
    p2_min=-500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=π / 4,
    sides_per_poly=4,
    derivs_per_sd=4, # 1
    derivs_per_fv=4, # 1
    n_obs=4,
    is_newsd=false
)
    n_x = 6
    n_u = 3
    n_xu = n_x + n_u

    n_ego = length(ego_polys)
    Ao = Symbolics.@variables(Ao[1:sides_per_poly, 1:2])[1] |> Symbolics.scalarize
    bo = Symbolics.@variables(bo[1:sides_per_poly])[1] |> Symbolics.scalarize
    centroido = Symbolics.@variables(centroido[1:2])[1] |> Symbolics.scalarize
    z = Symbolics.@variables(z[1:n_xu*T])[1] |> Symbolics.scalarize
    xt = Symbolics.@variables(xt[1:n_x])[1] |> Symbolics.scalarize
    λsd = Symbolics.@variables(λsd)
    α = Symbolics.@variables(α)
    β = Symbolics.@variables(β)
    x0 = Symbolics.@variables(x0[1:n_x])[1] |> Symbolics.scalarize

    sds = map(ego_polys) do P
        Ae = collect(P.A)
        be = P.b
        Ve = P.V
        centroide = sum(Ve) / length(Ve)
        g_col_single(xt, Ae, be, centroide, Ao, bo, centroido; is_newsd=is_newsd)
    end

    num_sd_cons = T * n_obs * n_ego
    num_sd_mults = T * n_obs * derivs_per_sd * n_ego

    if enable_fvals
        fvals = map(ego_polys) do P
            Ae = collect(P.A)
            be = P.b
            f_pack_single(xt, Ae, be, Q, q)
        end
        fkeys = map(fvals) do fv
            collect(keys(fv))
        end
        num_f_mults = derivs_per_fv * n_ego
    else
        num_f_mults = 0
    end

    cost_nom = f(z, T, Rf, Qf) #sum of interdimate cost for control and position
    cons_dyn = g_dyn(z, x0, T, dt)
    cons_env = g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)
    cons_nom = [cons_dyn; cons_env] #dynamic and environment constraints

    # coefficients of subderivatives of final cost
    if num_f_mults > 0
        α_f = Symbolics.@variables(α_f[1:num_f_mults])[1] |> Symbolics.scalarize
    else
        α_f = Vector{Symbolics.Num}()
    end
    # coefficients of subderivatives of sd (for every pair of ego and obstacle and every timestep)
    β_sd = Symbolics.@variables(β_sd[1:num_sd_mults])[1] |> Symbolics.scalarize
    # simplex constraints
    if num_f_mults > 0
        slacks = Symbolics.@variables(slacks[1:num_sd_cons+n_ego])[1] |> Symbolics.scalarize
    else
        slacks = Symbolics.@variables(slacks[1:num_sd_cons])[1] |> Symbolics.scalarize
    end
    # dynamic and environment constraints
    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize
    # collision avoidance constraints
    λ_col = Symbolics.@variables(λ_col[1:num_sd_cons])[1] |> Symbolics.scalarize

    θ = [z; α_f; β_sd; slacks; λ_nom; λ_col]

    lag = cost_nom - cons_nom' * λ_nom
    grad_lag = Symbolics.gradient(lag, z)

    simplex_cons = Num[]
    if num_f_mults > 0
        for i in 1:n_ego
            push!(simplex_cons, 1.0 - sum(α_f[(i-1)*derivs_per_fv+1:i*derivs_per_fv]))
        end
    end
    for i in 1:n_ego
        offset = (i - 1) * T * n_obs * derivs_per_sd
        for t in 1:T
            for j in 1:n_obs
                push!(simplex_cons, 1.0 - sum(β_sd[(1:derivs_per_sd).+(offset+(t-1)*n_obs*derivs_per_sd+(j-1)*derivs_per_sd)]))
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

    get_lag = Dict()
    get_sd = Dict()
    get_Jlag = Dict()
    get_Jsd = Dict()

    if enable_fvals
        get_gfv = Dict()
        get_Jfv = Dict()
    end

    @info "Generating symbolic solutions to sd calculation"
    for (i, sds_i) in enumerate(sds) # i is index of egos
        @showprogress for (k, sd) in sds_i # k is indices of active constraints, sd is corresponding sd
            lag = Symbolics.gradient(-sd * λsd[1] * β[1], xt; simplify=false) # λ is a representative of λ_col, β is a representative of β_sd
            get_lag[i, k] = Symbolics.build_function(lag, xt, Ao, bo, centroido, β, λsd; expression=Val(false))[2]
            get_sd[i, k] = Symbolics.build_function(β[1] * sd, xt, Ao, bo, centroido, β, λsd; expression=Val(false))

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
            get_Jlag[i, k] = (Jlag_rows, Jlag_cols, Symbolics.build_function(Jlag_vals, xt, Ao, bo, centroido, β, λsd; expression=Val(false))[2], zeros(length(Jlag_rows)))
            get_Jsd[i, k] = (Jsd_rows, Jsd_cols, Symbolics.build_function(Jsd_vals, xt, Ao, bo, centroido, β, λsd; expression=Val(false))[2], zeros(length(Jsd_rows)))
        end
    end
    if enable_fvals
        @info "Generating symbolic solutions to f calculation"
        for (i, fvals_i) in enumerate(fvals)
            @showprogress for (k, fv) in fvals_i
                gfv = Symbolics.gradient(fv * α[1], xt; simplify=false) # α is the coefficient of one subderivative
                get_gfv[i, k] = Symbolics.build_function(gfv, xt, α; expression=Val(false))[2]

                Jfv = Symbolics.sparsejacobian(gfv, [xt; α]; simplify=false)
                Jfv_rows, Jfv_cols, Jfv_vals = findnz(Jfv)
                get_Jfv[i, k] = (Jfv_rows, Jfv_cols, Symbolics.build_function(Jfv_vals, xt, α; expression=Val(false))[2], zeros(length(Jfv_rows)))
            end
        end
    end
    Aes = [deepcopy(P.A) for P in ego_polys]
    bes = [deepcopy(P.b) for P in ego_polys]
    centroides = [sum(deepcopy(P.V)) / length(deepcopy(P.V)) for P in ego_polys]

    col_offset = length(z) + length(α_f) + length(β_sd) + length(slacks) + length(λ_nom)
    β_offset = length(z) + length(α_f)
    lag_buf = zeros(n_x)

    function fill_F!(F, z, x0, polys, α_f, β_sd, λ_nom, λ_col)
        F .= 0.0
        get_Fnom!(F, z, x0, λ_nom, α_f, β_sd)

        for t in 1:T
            for i in 1:n_ego
                xt_inds = (t-1)*n_xu+1:(t-1)*n_xu+n_x
                @inbounds xt = z[xt_inds]
                for (e, Po) in enumerate(polys)
                    Ao = Po.A
                    bo = Po.b
                    centroido = sum(Po.V) / length(Po.V)
                    λ_ind = (i - 1) * T * n_obs + (t - 1) * n_obs + e
                    β_inds = (1:derivs_per_sd) .+ (((i - 1) * T * n_obs * derivs_per_sd) + (t - 1) * n_obs * derivs_per_sd + (e - 1) * derivs_per_sd)
                    @inbounds λte = λ_col[λ_ind]
                    @inbounds βte = β_sd[β_inds]

                    assignments = get_single_sd_ids(xt, Aes[i], bes[i], centroides[i], Ao, bo, centroido, derivs_per_sd, 0; is_newsd=is_newsd)
                    # directly delete indices after 3 may cause some problems
                    #if t == 20
                    @info "f assignments $assignments"
                    #end
                    #@infiltrate
                    #@infiltrate
                    for (ee, assignment) in enumerate(assignments)
                        if length(assignment) > 3
                            assignment = assignment[1:3]
                        end
                        get_lag[i, assignment](lag_buf, xt, Ao, bo, centroido, βte[ee], λte)
                        βsd = get_sd[i, assignment](xt, Ao, bo, centroido, βte[ee], λte)

                        F[xt_inds] .+= lag_buf
                        F[λ_ind+col_offset] += βsd
                    end
                end

                if enable_fvals && t == T
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

    function get_J_both!(JJ, z, x0, polys, α_f, β_sd, λ_nom, λ_col)
        JJ.nzval .= 1e-16
        get_Jnom_vals(Jnom_buf, z, x0, λ_nom, α_f, β_sd)
        JJ .+= sparse(Jnom_rows, Jnom_cols, Jnom_buf, n, n)

        for t in 1:T
            for i in 1:n_ego
                xt_inds = (t-1)*n_xu+1:(t-1)*n_xu+n_x
                @inbounds xt = z[xt_inds]
                for (e, P) in enumerate(polys)
                    Ao = collect(P.A)
                    bo = P.b
                    λ_ind = (i - 1) * T * n_obs + (t - 1) * n_obs + e
                    β_inds = (1:derivs_per_sd) .+ (((i - 1) * T * n_obs * derivs_per_sd) + (t - 1) * n_obs * derivs_per_sd + (e - 1) * derivs_per_sd)
                    @inbounds λte = λ_col[λ_ind]
                    @inbounds βte = β_sd[β_inds]
                    assignments = get_single_sd_ids(xt, Aes[i], bes[i], centroides[i], Ao, bo, centroido, derivs_per_sd, t; is_newsd=is_newsd)
                    @info "J assignments $assignments"
                    for (ee, assignment) in enumerate(assignments)
                        if length(assignment) > 3
                            assignment = assignment[1:3]
                        end
                        Jlag_rows, Jlag_cols, Jlag_vals, Jlag_buf = get_Jlag[i, assignment]
                        Jlag_vals(Jlag_buf, xt, Ao, bo, centroido, βte[ee], λte)
                        Jlag = sparse(Jlag_rows, Jlag_cols, Jlag_buf, 6, 8)
                        Jsd_rows, Jsd_cols, Jsd_vals, Jsd_buf = get_Jsd[i, assignment]
                        Jsd_vals(Jsd_buf, xt, Ao, bo, centroido, βte[ee], λte)
                        Jsd = sparse(Jsd_rows, Jsd_cols, Jsd_buf, 1, 7)

                        @inbounds JJ[xt_inds, xt_inds] .+= Jlag[1:6, 1:6]
                        @inbounds JJ[xt_inds, β_inds[ee]+β_offset] .+= Jlag[1:6, 7]
                        @inbounds JJ[xt_inds, λ_ind+col_offset] .+= Jlag[1:6, 8]
                        @inbounds JJ[λ_ind+col_offset, xt_inds] .+= Jsd[1:6]
                        @inbounds JJ[λ_ind+col_offset, β_inds[ee]+β_offset] += Jsd[7]
                    end
                end

                if enable_fvals && t == T
                    assignments_f = get_single_f_pack_ids(xt, Aes[i], bes[i], Q, q, derivs_per_fv, fkeys[i])
                    α_inds = (i-1)*derivs_per_fv+1:i*derivs_per_fv
                    @inbounds αi = α_f[α_inds]
                    for (ee, assignment) in enumerate(assignments_f)
                        Jfv_rows, Jfv_cols, Jfv_vals, Jfv_buf = get_Jfv[i, assignment]
                        Jfv_vals(Jfv_buf, xt, αi[ee])
                        Jfv = sparse(Jfv_rows, Jfv_cols, Jfv_buf, 6, 7)

                        @inbounds JJ[xt_inds, xt_inds] .+= Jfv[:, 1:length(xt_inds)]
                        @inbounds JJ[xt_inds, length(z)+α_inds[ee]] .+= Jfv[:, end]
                    end
                end
            end
        end
        nothing
    end

    J_example = sparse(Jnom_rows, Jnom_cols, ones(length(Jnom_cols)), n, n)

    for t in 1:T
        for i in 1:n_ego
            xt_inds = (t-1)*n_xu+1:(t-1)*n_xu+n_x
            for e in 1:n_obs
                λ_ind = (i - 1) * T * n_obs + (t - 1) * n_obs + e
                β_inds = (1:derivs_per_sd) .+ (((i - 1) * T * n_obs * derivs_per_sd) + (t - 1) * n_obs * derivs_per_sd + (e - 1) * derivs_per_sd)
                J_example[xt_inds, xt_inds] .= 1.0
                J_example[xt_inds, col_offset+λ_ind] .= 1.0
                J_example[col_offset+λ_ind, xt_inds] .= 1.0
                for l in β_inds
                    J_example[xt_inds, β_offset+l] .= 1.0
                    J_example[col_offset+λ_ind, β_offset+l] = 1.0
                end
            end
            if enable_fvals && t == T
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
    centroider = randn(sides_per_poly, 2)
    centroidor = randn(sides_per_poly, 2)
    lag_buf = zeros(6)
    λsdr = randn()
    αr = randn()
    βr = randn()
    for i in 1:n_ego
        @showprogress for k in collect(keys(sds[i]))
            get_lag[i, k](lag_buf, xtr, Aor, bor, centroidor, βr, λsdr)
            ss = get_sd[i, k](xtr, Aor, bor, centroidor, βr, λsdr)
            get_Jlag[i, k][3](get_Jlag[i, k][4], xtr, Aor, bor, centroidor, βr, λsdr)
            get_Jsd[i, k][3](get_Jsd[i, k][4], xtr, Aor, bor, centroidor, βr, λsdr)
            #ll = get_lag[k](xtr,Aer,ber,Aor,bor, λsdr)
            #ss = get_sd[k](xtr,Aer,ber,Aor,bor, λsdr)
            #J1 = get_Jlag[k](xtr,Aer,ber,Aor,bor,λsdr)
            #J2 = get_Jsd[k](xtr,Aer,ber,Aor,bor,λsdr)
        end
        if enable_fvals
            @showprogress for k in collect(keys(fvals[i]))
                gfv = get_gfv[i, k](lag_buf, xtr, αr)
                Jfv = get_Jfv[i, k][3](get_Jfv[i, k][4], xtr, αr)
            end
        end
    end

    return (; fill_F!,
        get_J_both!,
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
        n_obs,
        ego_polys,
        sides_per_poly,
        p1_max,
        p2_min)
end

function visualize_quick(x0, T, ego_polys, obs_polys; fig=Figure(), ax=Axis(fig[1, 1], aspect=DataAspect()), θ=[], is_displaying=true, is_newsd=false)
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

        # P = obs_polys[1]
        # Ao = P.A
        # bo = P.b
        # centroido = sum(P.V) / length(P.V)
        # LP_data = @lift(gen_LP_data($(xxts[i, t]), $(Aeb)[1], $(Aeb)[2], Ao, bo, centroido; is_newsd=is_newsd))
        # AA = @lift($(LP_data)[1])
        # bb = @lift($(LP_data)[2])
        # qq = @lift($(LP_data)[3])
        # ret = @lift(solve_qp(UseOSQPSolver(); A=sparse($(AA)), l=-$(bb), q=$(qq), polish=true, verbose=false))
        # sd = @lift($(ret).x[3])

        # @infiltrate
        # if !is_newsd
        #     be_inflated = @lift($(Aeb)[2] + $(sd) * ones(length($(Aeb)[2])))
        #     bo_inflated = @lift(bo + $(sd) * ones(length(bo)))            
        # else
        #     be_inflated = @lift($(Aeb)[2] + $(sd) * (($(Aeb)[1]) * $(xxts[i, t])[1:2] + ($(Aeb)[2])))
        #     bo_inflated = @lift(bo + $(sd) * (Ao * centroido + bo))
        # end
        # self_poly_inflated = @lift(ConvexPolygon2D($(Aeb)[1], $(be_inflated)))
        # plot!(ax, self_poly_inflated; color=:yellow, linewidth=3)

        # obstacle_inflated = @lift(ConvexPolygon2D(Ao, $(bo_inflated)))
        # plot!(ax, obstacle_inflated; color=:yellow, linewidth=3)
    end

    #colors = [:red, :orange, :yellow, :green]
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

function solve_quick(prob, x0, obs_polys; θ0=nothing, is_displaying=true, sleep_duration=0.0, is_newsd=false)
    (; fill_F!, get_J_both!, J_example, ego_polys, l, u, T, n_z, n_α, n_β, n_s, n_nom, n_col, n_obs, sides_per_poly, p1_max, p2_min) = prob

    @assert length(obs_polys) == n_obs

    n = length(l)
    @assert n == n_z + n_α + n_β + n_s + n_nom + n_col

    if is_displaying
        (fig, update_fig) = visualize_quick(x0, T, ego_polys, obs_polys; is_newsd)
    end

    if isnothing(θ0)
        θ0 = zeros(n)
        derivs_per_fv = Int(n_α / length(ego_polys))
        derivs_per_sd = Int(n_β / (T * length(ego_polys) * n_obs))
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
        @inbounds z = θ[1:n_z]
        @inbounds α_f = θ[n_z+1:n_z+n_α]
        @inbounds β_sd = θ[n_z+n_α+1:n_z+n_α+n_β]
        @inbounds λ_nom = θ[n_z+n_α+n_β+n_s+1:n_z+n_α+n_β+n_s+n_nom]
        @inbounds λ_col = θ[n_z+n_α+n_β+n_s+n_nom+1:n_z+n_α+n_β+n_s+n_nom+n_col]
        fill_F!(result, z, x0, obs_polys, α_f, β_sd, λ_nom, λ_col)

        if is_displaying
            update_fig(θ)
            if sleep_duration > 0
                sleep(sleep_duration)
            end
        end

        Cint(0)
    end
    function F2(n, θ, result)
        result .= 0.0
        @inbounds z = θ[1:n_z]
        @inbounds α_f = θ[n_z+1:n_z+n_α]
        @inbounds β_sd = θ[n_z+n_α+1:n_z+n_α+n_β]
        @inbounds λ_nom = θ[n_z+n_α+n_β+n_s+1:n_z+n_α+n_β+n_s+n_nom]
        @inbounds λ_col = θ[n_z+n_α+n_β+n_s+n_nom+1:n_z+n_α+n_β+n_s+n_nom+n_col]
        fill_F!(result, z, x0, obs_polys, α_f, β_sd, λ_nom, λ_col)

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
        get_J_both!(JJ, z, x0, obs_polys, α_f, β_sd, λ_nom, λ_col)

        col .= J_col
        len .= J_len
        row .= J_row
        data .= JJ.nzval

        @info "filling base F"
        buf = zeros(n)
        buf2 = zeros(n)
        F(n, θ, buf)
        @info "filling perturbed F"
        #Jrows, Jcols, _ = findnz(J_example)

        #Main.@infiltrate
        #Jnum = sparse(Jrows, Jcols, Jbuf)
        Jnum2 = spzeros(n, n)
        #@info "Testing Jacobian accuracy numerically"
        @showprogress for ni in 1:n
            wi = copy(θ)
            wi[ni] += 1e-3
            F2(n, wi, buf2)
            Jnum2[:, ni] = sparse((buf2 - buf) ./ 1e-3)
        end
        @info "Jacobian error is $(norm(Jnum2-JJ))"
        if norm(Jnum2 - JJ) > 1e-1
            @infiltrate
        end


        #if is_displaying
        #    update_fig(θ)
        #    if sleep_duration > 0
        #        sleep(sleep_duration)
        #    end
        #end

        Cint(0)
    end

    # force compilation
    buf = zeros(n)
    Jbuf = zeros(nnz_total)
    w = randn(length(θ0))
    F(n, w, buf)
    J(n, nnz_total, w, zero(J_col), zero(J_len), zero(J_row), Jbuf)

    #buf2 = zeros(n)
    #Jrows, Jcols, _ = findnz(J_example)
    #Main.@infiltrate
    #Jnum = sparse(Jrows, Jcols, Jbuf)
    #Jnum2 = spzeros(n, n)
    #@info "Testing Jacobian accuracy numerically"
    #@showprogress for ni in 1:n
    #    wi = copy(w)
    #    wi[ni] += 1e-5
    #    F(n, wi, buf2)
    #    Jnum2[:,ni] = sparse((buf2-buf) ./ 1e-5)
    #end
    #@info "Jacobian error is $(norm(Jnum2-Jnum))"

    PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")
    status, θ, info = PATHSolver.solve_mcp(
        F,
        J,
        l,
        u,
        θ0;
        silent=false,
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

    if is_displaying
        display(fig)
    end

    @inbounds z = @view(θ[1:n_z])

    (; status, info, θ, z, fres)
end




