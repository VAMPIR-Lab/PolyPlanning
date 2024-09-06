
function setup_nonsmooth_3d_0step(
    ego_polys,
    obs_polys;
    T=2,
    Q_cost=1e-3 * I(3), # penality for distance
    p1_max=500.0,
    p2_max=500.0,
    p3_max=500.0,
    n_sd_slots=2
)

    # problem dimensions
    n_x = 6
    n_u = 0
    n_xu = n_x + n_u
    n_z = T * n_xu
    n_ego = length(ego_polys)
    n_obs = length(obs_polys)
    # n_side_ego = length(ego_polys[1].b)
    # n_side_obs = length(obs_polys[1].b)
    n_side_ego = maximum([length(i.b) for i in ego_polys]) # just for the buffer
    n_side_obs = maximum([length(i.b) for i in obs_polys]) # just for the buffer
    n_env_cons = T * n_xu
    # combin_2_from_n = n::Int -> n * (n - 1) ÷ 2
    # n_sds = (combin_2_from_n(n_side_ego) * n_side_obs + n_side_ego * combin_2_from_n(n_side_obs))
    n_sd_cons = T * n_ego * n_obs * n_sd_slots

    # assume number of sides are the same for all ego polys, obs polys
    # @assert all([length(ego_polys[i].b) for i in 1:n_ego] .== n_side_ego)
    # @assert all([length(obs_polys[i].b) for i in 1:n_obs] .== n_side_obs)

    # θ indexing
    z_s2i, env_cons_s2i, sd_cons_s2i = get_sub2idxs((n_xu, T), (n_env_cons), (n_sd_slots, n_obs, n_ego, T))

    z = Symbolics.@variables(z[1:n_z])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:n_x])[1] |> Symbolics.scalarize
    xt = Symbolics.@variables(xt[1:n_x])[1] |> Symbolics.scalarize

    cost = f_3d_0step(z, T, Q_cost)
    env_cons = g_env_3d_0step(z, T, p1_max, p2_max, p3_max)

    # check indexing consistency
    @assert length(z_s2i) == n_z
    @assert length(env_cons_s2i) == length(env_cons)
    @assert length(sd_cons_s2i) == n_sd_cons

    # sds, intercepts, AAs, bbs, etc to fill F and J
    get_AA = Dict()
    get_bb = Dict()
    Itr_dict = Dict()

    # functions vary according to different egos
    get_grad_A_ego_wrt_xt = Dict()
    get_grad_b_ego_wrt_xt = Dict()
    get_Hessian_A_ego_wrt_xt = Dict()
    get_Hessian_b_ego_wrt_xt = Dict()
    # we solve sds symbolically for given ego and obs at problem creation, 
    # TODO could be done in a more flexible by abstracting problem parameters (obs_polys) and filling them in fill_F, fill_J instead as we did before
    λsd = Symbolics.@variables(λsd)[1]

    # rotation matrix and translation 
    R = R_from_mrp(xt[4:6])
    l = xt[1:3]

    # helper functions to calculate Jacobian
    aibi_wrt_xt_functions = get_aibi_wrt_xt_fun(R, l, xt[1:6])
    get_grad_sd_wrt_Ab, get_Hessian_sd_wrt_Ab, get_sol_A_b = get_sd_wrt_Ab_fun()

    #@info "Generating symbolic sds, intercepts, AAs, bbs"
    #p = Progress(n_sds * n_ego * n_obs, dt=1.0)
    for (i, Pe) in enumerate(ego_polys)
        A_ego = collect(Pe.A)
        b_ego = Pe.b
        V_ego = Pe.V
        centr_ego = Pe.c
        get_grad_A_ego_wrt_xt[i], get_grad_b_ego_wrt_xt[i], get_Hessian_A_ego_wrt_xt[i], get_Hessian_b_ego_wrt_xt[i] = get_Ab_ego_wrt_xt_fun(xt[1:6], A_ego, b_ego, aibi_wrt_xt_functions)

        for (j, Po) in enumerate(obs_polys)
            A_obs = collect(Po.A)
            b_obs = Po.b
            V_obs = Po.V
            centr_obs = Po.c
            
            AA = [A_ego*R' A_ego*centr_ego+b_ego
                A_obs    A_obs*centr_obs+b_obs]
            bb = [b_ego-A_ego*R'*l
                b_obs]

            Itr_dict[i, j] = get_possible_assignments_3d(A_ego, b_ego, V_ego, A_obs, b_obs, V_obs)
            get_AA[i, j] = Symbolics.build_function(AA', xt; expression=Val(false))[2]
            get_bb[i, j] = Symbolics.build_function(bb, xt; expression=Val(false))[2]

        end
    end
    # different ego-obs pair may have different numbers of possible assignments
    n_sds = []
    for (i, dic) in enumerate(Itr_dict)
        push!(n_sds, length(dic[2]))
    end
    n_sds = maximum(n_sds)

    nom_cons = env_cons
    λ_nom = Symbolics.@variables(λ_nom[1:length(nom_cons)])[1] |> Symbolics.scalarize
    λ_sd = Symbolics.@variables(λ_col[1:n_sd_cons])[1] |> Symbolics.scalarize

    # F = [F_nom; F_sd (to be filled later)]
    # J = dJ/dθ
    # we also need sds/dx

    # F and J nominal (before filling)
    θ = [z; λ_nom; λ_sd]
    lag = cost - nom_cons' * λ_nom #- λ_sd' * nom_sd (to be filled later)
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; nom_cons; zeros(Num, n_sd_cons)]
    J_nom = Symbolics.sparsejacobian(F_nom, θ)

    n = length(F_nom)
    l = zeros(n)
    u = zeros(n)

    # gradient = 0
    l[z_s2i[:]] .= -Inf
    u[z_s2i[:]] .= Inf

    # environmental constraints
    l[env_cons_s2i[:]] .= 0.0
    u[env_cons_s2i[:]] .= Inf

    # sd constraints
    l[sd_cons_s2i[:]] .= 0.0
    u[sd_cons_s2i[:]] .= Inf

    # get nominal F and J (later to be filled by fill_F! and fill_J!)
    get_Fnom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false))[2]
# @infiltrate
    (Jnom_rows, Jnom_cols, Jnom_vals) = findnz(J_nom)
    get_Jnom_vals! = Symbolics.build_function(Jnom_vals, z, x0, λ_nom; expression=Val(false))[2]

    # sort sds and intercepts for given xt
    sol_xt_buffer_full = zeros(4, n_sds)
    AA_buffer_full = zeros(4, n_side_ego + n_side_obs)
    bb_buffer_full = zeros(n_side_ego + n_side_obs)

    function get_sorted_sds_3d(i, m1, j, m2, xt; tol=1e-4, local_factor=1.5)

        # if the size of buffer is larger than what the function returns, it is filled by column, so here AA is transposed
        get_AA[i, j](AA_buffer_full, xt)
        get_bb[i, j](bb_buffer_full, xt)
        AA_buffer = AA_buffer_full'[1:m1+m2, :]
        bb_buffer = bb_buffer_full[1:m1+m2]

        assignments = Itr_dict[i, j]
        n_ass = length(assignments)
        get_sol_xt!(sol_xt_buffer_full, AA_buffer, bb_buffer, get_sol_A_b, assignments)
        sol_xt_buffer = @view sol_xt_buffer_full'[1:n_ass, :]
        
        sds_buffer = @view sol_xt_buffer[:, end]
        intercept_buffer = @view sol_xt_buffer[:, 1:end-1]

        # tol = 1e-4
        # sd must be greater than -1 to be valid
        valid_mask = sds_buffer .>= -1.0 - tol

        # sd and intercept must correspond to a feasible vertex to be valid
        zz_check = hcat(intercept_buffer[valid_mask, :], sds_buffer[valid_mask])

        valid_mask[valid_mask] = map(eachrow(zz_check)) do row
            all(AA_buffer * row + bb_buffer .>= -tol)
        end

        sorted_sds_inds = sds_buffer[valid_mask] |> sortperm
        sorted_sds = sds_buffer[valid_mask][sorted_sds_inds]
        sorted_ass = assignments[valid_mask][sorted_sds_inds]

        # need to consider smarter way to filter out distant sds
        # local_factor = 1.5 # regard sds which are less than sd*local_factor as potential true sds
        local_sd_mask = (sorted_sds .+ 1) .<= (sorted_sds[1] + 1) * local_factor
        sorted_sds = sorted_sds[local_sd_mask]
        sorted_ass = sorted_ass[local_sd_mask]
        (sorted_sds, sorted_ass, AA_buffer, bb_buffer)
    end

    # compute the most consistent mapping between sd slot memory and current assignments
    sd_slot_mem = [Vector{Int64}() for _ in 1:n_sd_slots]

    function compute_ass_ind_map(sorted_ass, n_ass)
        ind_dict = Dict()
        ind_map = collect(1:n_ass)

        # identify which assignments exist in sd_slot_mem
        #for i in 1:n_ass
        #    for (k, mem) in enumerate(sd_slot_mem)
        #        if sorted_ass[i] == mem
        #            ind_dict[i] = k
        #        end
        #    end
        #end

        ## place remaining keys arbitrarily
        #for i in 1:n_ass
        #    if !haskey(ind_dict, i)
        #        k = 1
        #        while k ∈ values(ind_dict)
        #            k += 1
        #        end

        #        ind_dict[i] = k
        #    end
        #end

        #for (k, i) in ind_dict
        #    ind_map[k] = i
        #end

        #@infiltrate any(ind_map .!= collect(1:n_ass))
        ind_map = collect(1:n_ass) # smart indexing disabled
        ind_map
    end

    # fill_F!
    λ_nom_s2i = [env_cons_s2i...]

    # do not use fill(), because every element points to the same vector zeros(6)
    grad_A_ego_wrt_xt_buffer = [zeros(6) for _ in 1:n_side_ego, _ in 1:3]
    grad_b_ego_wrt_xt_buffer = [zeros(6) for _ in 1:n_side_ego]
    Hessian_A_ego_wrt_xt_buffer = [zeros(6, 6) for _ in 1:n_side_ego, _ in 1:3]
    Hessian_b_ego_wrt_xt_buffer = [zeros(6, 6) for _ in 1:n_side_ego]

    grad_A_wrt_xt_buffer = [zeros(6) for _ in 1:4, _ in 1:4]
    grad_b_wrt_xt_buffer = [zeros(6) for _ in 1:4]
    Hessian_A_wrt_xt_buffer = [zeros(6, 6) for _ in 1:4, _ in 1:4]
    Hessian_b_wrt_xt_buffer = [zeros(6, 6) for _ in 1:4]

    Ab_buffer = zeros(20)
    grad_Ab_wrt_xt_buffer = [zeros(6) for _ in 1:20]
    Hessian_Ab_wrt_xt_buffer = [zeros(6, 6) for _ in 1:20]
    Jacobian_Ab_wrt_xt_buffer = zeros(20, 6)
    HJ_buffer = zeros(20, 6)

    grad_sd_wrt_Ab_buffer = zeros(20)
    Hessian_sd_wrt_Ab_buffer = zeros(20, 20)

    grad_sd_wrt_xt_buffer = zeros(6)
    Hessian_sd_wrt_xt_buffer = zeros(6, 6)

    sd_lag_buf = zeros(n_x)
    Jsdlag_buf = zeros(n_x, n_x+1)
    Jsd_buf = zeros(n_x)

    function fill_F!(F, θ, x0)
        # TODO obs_polys as parameters
        F .= 0.0 # clear
        @inbounds z = θ[z_s2i[:]]
        @inbounds λ_nom = θ[λ_nom_s2i[:]]
        get_Fnom!(F, z, x0, λ_nom)

        for t in 1:T
            xt_ind = z_s2i[1:n_x, t]
            @inbounds xt = z[xt_ind]
            # println("F ", t)
            for (i, Pe) in enumerate(ego_polys)
                m1 = length(Pe.b)

                get_grad_A_ego_wrt_xt[i](grad_A_ego_wrt_xt_buffer, xt[1:6])
                get_grad_b_ego_wrt_xt[i](grad_b_ego_wrt_xt_buffer, xt[1:6])
                # get_J_A_ego_wrt_xt_fun[i](J_A_ego_wrt_xt_buffer, xt[1:6])
                # get_J_b_ego_wrt_xt_fun[i](J_b_ego_wrt_xt_buffer, xt[1:6])
                for (j, Po) in enumerate(obs_polys)
                    (sorted_sds, sorted_ass, AA, bb) = get_sorted_sds_3d(i, length(Pe.b), j, length(Po.b), xt)
                    n_ass = min(length(sorted_ass), n_sd_slots)
                    k_map = compute_ass_ind_map(sorted_ass, n_ass)


                    for slot_i in 1:n_sd_slots
                        if slot_i < n_ass
                            # smart allocation for existing assignments
                            sd_rank = slot_i
                            k = k_map[sd_rank]
                        else
                            # copy last assignment for others
                            sd_rank = n_ass
                            k = slot_i
                        end
                        ass = sorted_ass[sd_rank]

                        sd_ind = sd_cons_s2i[k, j, i, t]
                        @inbounds λsd = θ[sd_ind]
                        
                        get_grad_A_b_wrt_xt!(grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer, ass, grad_A_ego_wrt_xt_buffer, grad_b_ego_wrt_xt_buffer, m1)
                        # get_Hessian_A_b_wrt_xt!(Hessian_A_wrt_xt_buffer, Hessian_b_wrt_xt_buffer, ass, Hessian_A_ego_wrt_xt_buffer, Hessian_b_ego_wrt_xt_buffer, m1)

                        # Ab = vcat(vec(AA[ass,:]), bb[ass])
                        get_Ab!(Ab_buffer, AA[ass,:], bb[ass])
                        get_grad_sd_wrt_Ab(grad_sd_wrt_Ab_buffer, Ab_buffer)
                        # get_Hessian_sd_wrt_Ab(Hessian_sd_wrt_Ab_buffer, Ab_buffer)

                        get_grad_Ab_wrt_xt(grad_Ab_wrt_xt_buffer, grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer) # rearrange, get a vector of gradient vectorss
                        # get_Hessian_Ab_wrt_xt(Hessian_Ab_wrt_xt_buffer, Hessian_A_wrt_xt_buffer, Hessian_b_wrt_xt_buffer) # rearrange, get a vector of Hessian matrices
                        # get_Jacobian_Ab_wrt_xt(Jacobian_Ab_wrt_xt_buffer, grad_Ab_wrt_xt_buffer) # convert a vector of vectors to a matrix
                        # get_HJ(HJ_buffer, Hessian_sd_wrt_Ab_buffer, Jacobian_Ab_wrt_xt_buffer) # HJ = H * J

                        get_grad_sd_wrt_xt!(grad_sd_wrt_xt_buffer, grad_sd_wrt_Ab_buffer, grad_Ab_wrt_xt_buffer)
                        # get_Hessian_sd_wrt_xt!(Hessian_sd_wrt_xt_buffer, grad_sd_wrt_Ab_buffer, Jacobian_Ab_wrt_xt_buffer, HJ_buffer, Hessian_Ab_wrt_xt_buffer)

                        sd_lag_buf .= -λsd * grad_sd_wrt_xt_buffer

                        @inbounds F[xt_ind] += sd_lag_buf
                        @inbounds F[sd_ind] += sorted_sds[sd_rank]
                    end

                    # not updating sd slot memory in F is less robust
                    #sd_slot_mem[1:n_ass] = sorted_ass[1:n_ass]
                end
            end
        end
        nothing
    end

    # fill_J!
    Jnom_buf = zeros(length(Jnom_vals))

    function fill_J_vals!(J_vals, θ, x0)
        ### check Jacobian numerically
        #buf = zeros(n)
        #fill_F!(buf, θ, x0)
        #buf2 = zeros(n)
        #Jnum = spzeros(n, n)
        #for ni in 1:n
        #    wi = deepcopy(θ)
        #    wi[ni] += 1e-5
        #    fill_F!(buf2, wi, x0)
        #    #@infiltrate ni == 49
        #    Jnum[:, ni] = sparse((buf2 - buf) ./ 1e-5)
        #end
        ####

        # TODO obs_polys as parameters
        J_vals.nzval .= 1e-16 # clear
        @inbounds z = θ[z_s2i[:]]
        @inbounds λ_nom = θ[λ_nom_s2i[:]]

        get_Jnom_vals!(Jnom_buf, z, x0, λ_nom)
        J_vals .+= sparse(Jnom_rows, Jnom_cols, Jnom_buf, n, n)

        for t in 1:T
            xt_ind = z_s2i[1:n_x, t]
            @inbounds xt = z[xt_ind]
            # println("J ", t)
            for (i, Pe) in enumerate(ego_polys)
                m1 = length(Pe.b)
# @infiltrate
                get_grad_A_ego_wrt_xt[i](grad_A_ego_wrt_xt_buffer, xt[1:6])
                get_grad_b_ego_wrt_xt[i](grad_b_ego_wrt_xt_buffer, xt[1:6])
                get_Hessian_A_ego_wrt_xt[i](Hessian_A_ego_wrt_xt_buffer, xt[1:6])
                get_Hessian_b_ego_wrt_xt[i](Hessian_b_ego_wrt_xt_buffer, xt[1:6])
# @infiltrate
                for (j, Po) in enumerate(obs_polys)
                    (sorted_sds, sorted_ass, AA, bb) = get_sorted_sds_3d(i, length(Pe.b), j, length(Po.b), xt)

                    n_ass = min(length(sorted_ass), n_sd_slots)
                    k_map = compute_ass_ind_map(sorted_ass, n_ass)

                    for slot_i in 1:n_sd_slots
                        if slot_i < n_ass
                            # smart allocation for existing assignments
                            sd_rank = slot_i
                            k = k_map[sd_rank]
                        else
                            # copy last assignment for others
                            sd_rank = n_ass
                            k = slot_i
                        end
                        ass = sorted_ass[sd_rank]

                        sd_ind = sd_cons_s2i[k, j, i, t]
                        @inbounds λsd = θ[sd_ind]

                        get_grad_A_b_wrt_xt!(grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer, ass, grad_A_ego_wrt_xt_buffer, grad_b_ego_wrt_xt_buffer, m1)
                        get_Hessian_A_b_wrt_xt!(Hessian_A_wrt_xt_buffer, Hessian_b_wrt_xt_buffer, ass, Hessian_A_ego_wrt_xt_buffer, Hessian_b_ego_wrt_xt_buffer, m1)

                        # Ab = vcat(vec(AA[ass,:]), bb[ass])
                        get_Ab!(Ab_buffer, AA[ass,:], bb[ass])
                        get_grad_sd_wrt_Ab(grad_sd_wrt_Ab_buffer, Ab_buffer)
                        get_Hessian_sd_wrt_Ab(Hessian_sd_wrt_Ab_buffer, Ab_buffer)

                        get_grad_Ab_wrt_xt(grad_Ab_wrt_xt_buffer, grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer) # rearrange, get a vector of gradient vectorss
                        get_Hessian_Ab_wrt_xt(Hessian_Ab_wrt_xt_buffer, Hessian_A_wrt_xt_buffer, Hessian_b_wrt_xt_buffer) # rearrange, get a vector of Hessian matrices
                        get_Jacobian_Ab_wrt_xt(Jacobian_Ab_wrt_xt_buffer, grad_Ab_wrt_xt_buffer) # convert a vector of vectors to a matrix
                        get_HJ(HJ_buffer, Hessian_sd_wrt_Ab_buffer, Jacobian_Ab_wrt_xt_buffer) # HJ = H * J

                        get_grad_sd_wrt_xt!(grad_sd_wrt_xt_buffer, grad_sd_wrt_Ab_buffer, grad_Ab_wrt_xt_buffer)
                        get_Hessian_sd_wrt_xt!(Hessian_sd_wrt_xt_buffer, grad_sd_wrt_Ab_buffer, Jacobian_Ab_wrt_xt_buffer, HJ_buffer, Hessian_Ab_wrt_xt_buffer)

                        Jsd_buf .= grad_sd_wrt_xt_buffer
                        Jsdlag_buf .= hcat(-λsd*Hessian_sd_wrt_xt_buffer, -grad_sd_wrt_xt_buffer)


                        # the same as code block above
                        # get_Jsdlag!(Jsdlag_buf, xt, λsd, i, ass, AA, bb, m1)
                        # get_Jsd!(Jsd_buf, xt, i, ass, AA, bb, m1)

                        @inbounds J_vals[xt_ind, xt_ind] += Jsdlag_buf[1:n_x, 1:n_x]
                        @inbounds J_vals[xt_ind, sd_ind] += Jsdlag_buf[1:n_x, n_x+1]
                        @inbounds J_vals[sd_ind, xt_ind] += Jsd_buf
                        # if maximum(J_vals) > 1e8
                        #     @infiltrate
                        #     # plot_xt(xt, Pe.A, Pe.b, Po.A, Po.b)
                        # end
                    end

                    # update sd slot memory
                    sd_slot_mem[1:n_ass] = sorted_ass[1:n_ass]
                end
            end
        end
        # check if Jacobian has exploded
        # @infiltrate
        # max_val, max_ind = findmax(J_vals)
        # min_val, min_ind = findmin(J_vals)
        # println("max", Tuple(max_ind), "=", round(max_val; sigdigits=3), " min", Tuple(min_ind), "=", round(min_val; sigdigits=3))
        # println("J[30,494]", J_vals[30,494])
        nothing
    end

    J_example = sparse(Jnom_rows, Jnom_cols, ones(length(Jnom_cols)), n, n)

    for t in 1:T
        xt_ind = z_s2i[1:n_x, t]
        @inbounds xt = z[xt_ind]

        for (i, Pe) in enumerate(ego_polys)
            for (j, Po) in enumerate(obs_polys)
                sd_ind = sd_cons_s2i[:, j, i, t]
                @inbounds J_example[xt_ind, xt_ind] .+= 1.0
                @inbounds J_example[xt_ind, sd_ind] .+= 1.0
                @inbounds J_example[sd_ind, xt_ind] .+= 1.0
            end
        end
    end

    # forced compilation and J example
    # TODO verify this is okay
    #@info "Forcing compilation by computing J example"
    #J_example = sparse(Jnom_rows, Jnom_cols, ones(length(Jnom_cols)), n, n)
    #fill_J_vals!(J_example, rand(length(θ)), rand(length(x0)))

    param = (;
        ego_polys,
        obs_polys,
        T,
        Q_cost,
        p1_max,
        p2_max,
        p3_max,
        n_sd_slots,
        z_s2i,
        env_cons_s2i,
        sd_cons_s2i
    )

    (;
        fill_F!,
        fill_J_vals!,
        J_example,
        l,
        u,
        param
    )
end

# this needs to be updated for 3d
function visualize_nonsmooth_3d_0step(x0, T, ego_polys, obs_polys; fig=Figure(), ax3=LScene(fig[1, 1], scenekw=(camera=cam3d!, show_axis=true)), θ=[], is_displaying=true)
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    xxts = Dict()
    

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

function solve_nonsmooth_3d_0step(prob, x0; θ0=nothing, is_displaying=true, sleep_duration=0.0)
    param = prob.param

    n_x = 6
    n_u = 0
    n_xu = n_x + n_u
    n = length(prob.l)

    if is_displaying
        (fig, update_fig) = visualize_nonsmooth_3d_0step(x0, param.T, param.ego_polys, param.obs_polys)
    end

    # initialize
    if isnothing(θ0)
        θ0 = zeros(n)
        for t in 1:param.T
            θ0[(t-1)*n_xu+1:(t-1)*n_xu+n_x] = x0
        end
    end

    # F
    function F(n, θ, FF)
        FF .= 0.0
        prob.fill_F!(FF, θ, x0)

        if is_displaying
            update_fig(θ)
            if sleep_duration > 0
                sleep(sleep_duration)
            end
        end
        Cint(0)
    end

    # J
    J_vals = prob.J_example
    J_col = J_vals.colptr[1:end-1]
    J_len = diff(J_vals.colptr)
    J_row = J_vals.rowval
    nnz_total = length(J_vals.nzval)

    function J(n, nnz, θ, col, len, row, data)
        @assert nnz == nnz_total
        data .= 0.0

        prob.fill_J_vals!(J_vals, θ, x0)
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
    # t=time()
    J(n, nnz_total, w, zero(J_col), zero(J_len), zero(J_row), Jbuf)
    # println("update J takes ", time()-t, " seconds")
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
    #    if norm(buf2 - buf)>1e-3
    #     @infiltrate
    #    end
    # end
    # @info "Jacobian error is $(norm(Jnum2-Jnum))"
    PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")
    status, θ, info = PATHSolver.solve_mcp(
        F,
        J,
        prob.l,
        prob.u,
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

    f_res = zeros(n)
    F(n, θ, f_res)

    if is_displaying
        display(fig)
    end

    @inbounds z = @view(θ[param.z_s2i[:]])

    (; status, info, θ, z, f_res)
end