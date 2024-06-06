function get_lagrangian(xt, Ae, be, Ao, bo, p, sd, λsd)
    # min c'y
    # s.t. h - G y in K
    # KKT:
    # c + G'λ = 0
    # h - G y = 0
    # λ in K*
    # (h - G y) ∘ λ = 0

    # for polytopes:
    # min sd
    # s.t. 
    # A r - [A -b] [p, sd] ≥ R+
    # α ≥ 0
    # c' x = s + λ (cons)
    # p = 2 dim sym

    Aex, bex = shift_to(Ae, be, xt)
    c = [0; 0; 1.0]
    G = [Aex -bex; Ao -bo]
    h = [Aex * xt[1:2]; Ao * xt[1:2]]

    c' * [p; sd] + λsd' * (G * [p; sd] - h)
end

function setup_differentiablecol(ego_polys;
    T=50,
    dt=0.2,
    Rf=1e-3 * I(3),
    Qf=1e-3 * I(2),
    p1_max=500.0,
    p2_min=-500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=π / 4,
    sides_per_poly=4,
    n_obs=4
)
    n_x = 6
    n_u = 3
    n_xu = n_x + n_u

    n_ego = length(ego_polys)
    sides_per_ego = length(ego_polys[1].b)
    # assert all sides of polys are the same for all i
    @assert all([length(ego_polys[i].b) for i in 1:n_ego] .== sides_per_ego)

    Ao = Symbolics.@variables(Ao[1:sides_per_poly, 1:2])[1] |> Symbolics.scalarize
    bo = Symbolics.@variables(bo[1:sides_per_poly])[1] |> Symbolics.scalarize

    z = Symbolics.@variables(z[1:n_xu*T])[1] |> Symbolics.scalarize
    xt = Symbolics.@variables(xt[1:n_x])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:n_x])[1] |> Symbolics.scalarize

    # intersection point, sign distance, corresponding duals
    p = Symbolics.@variables(p[1:2])[1] |> Symbolics.scalarize
    sd = Symbolics.@variables(sd)
    λsd = Symbolics.@variables(λsd[1:sides_per_poly+sides_per_ego])[1] |> Symbolics.scalarize

    col_lags = map(ego_polys) do P
        Ae = collect(P.A)
        be = P.b
        get_lagrangian(xt, Ae, be, Ao, bo, p, sd, λsd)
    end

    cost_nom = f(z, T, Rf, Qf)
    cons_dyn = g_dyn(z, x0, T, dt)
    cons_env = g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)
    cons_nom = [cons_dyn; cons_env]

    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize
    θ = [z; λ_nom]

    lag = cost_nom - cons_nom' * λ_nom
    grad_lag = Symbolics.gradient(lag, z)

    F_nom = [grad_lag; cons_nom]

    l = [fill(-Inf, length(grad_lag)); fill(-Inf, length(cons_dyn)); zeros(length(cons_env))]
    u = [fill(+Inf, length(grad_lag)); fill(Inf, length(cons_dyn)); fill(Inf, length(cons_env))]
    n = length(l)

    Jnom = Symbolics.sparsejacobian(F_nom, θ)

    (Jnom_rows, Jnom_cols, Jnom_vals) = findnz(Jnom)
    Jnom_buf = zeros(length(Jnom_vals))
    get_Fnom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false))[2]
    get_Jnom_vals = Symbolics.build_function(Jnom_vals, z, x0, λ_nom; expression=Val(false))[2]

    get_dlag = Dict()
    get_Jlag = Dict()

    for (i, lag) in enumerate(col_lags)
        dlag = Symbolics.gradient(lag, xt; simplify=false)
        get_dlag[i] = Symbolics.build_function(dlag, xt, Ao, bo, λsd, p, sd; expression=Val(false))[2]
        Jlag = Symbolics.sparsejacobian(dlag, xt; simplify=false)

        Jlag_rows, Jlag_cols, Jlag_vals = findnz(Jlag)
        if length(Jlag_vals) == 0
            Jlag_rows = [1]
            Jlag_cols = [1]
            Jlag_vals = Num[0.0]
        end

        get_Jlag[i] = (Jlag_rows, Jlag_cols, Symbolics.build_function(Jlag_vals, xt, Ao, bo, λsd, p, sd; expression=Val(false))[2], zeros(length(Jlag_rows)))
    end

    #col_offset = length(z) + length(λ_nom)
    #lag_buf = zeros(6)
    Aes = [deepcopy(P.A) for P in ego_polys]
    bes = [deepcopy(P.b) for P in ego_polys]
    lag_buf = zeros(n_x)


    function fill_F!(F, z, x0, polys, λ_nom)
        F .= 0.0
        get_Fnom!(F, z, x0, λ_nom)
        for i in 1:n_ego
            for t in 1:T
                xt_inds = (t-1)*n_xu+1:(t-1)*n_xu+n_x
                @inbounds xt = z[xt_inds]
                for (e, P) in enumerate(polys)
                    Ao = P.A
                    bo = P.b

                    # osqp solver -> sd, p, λsd
                    Aex, bex = shift_to(Aes[i], bes[i], xt)
                    qq = [0; 0; 1.0]
                    AA = [Aex -bex; Ao -bo]
                    bb = [Aex * xt[1:2]; Ao * xt[1:2]]
                    # -AA * [p; sd] ≥ -bb
                    ret = solve_qp(UseOSQPSolver(); A=-sparse(AA), l=-bb, q=qq, polish=true, verbose=false)

                    p = ret.x[1:2]
                    sd = ret.x[3]
                    λsd = ret.y

                    get_dlag[i](lag_buf, xt, Ao, bo, λsd, p, sd)
                    F[xt_inds] .+= lag_buf
                end
            end
        end
        nothing
    end

    function get_J_both!(JJ, z, x0, polys, λ_nom)
        JJ = sparse(Jnom_rows, Jnom_cols, ones(length(Jnom_cols)), n, n)
        #JJ .= 0.0
        #JJ.nzval .= 1e-16
        #get_Jnom_vals(Jnom_buf, z, x0, λ_nom, α_f, β_sd)
        #JJ .+= sparse(Jnom_rows, Jnom_cols, Jnom_buf, n, n)

        for i in 1:n_ego
            for t in 1:T
                xt_inds = (t-1)*n_xu+1:(t-1)*n_xu+n_x
                @inbounds xt = z[xt_inds]
                for (e, P) in enumerate(polys)
                    Ao = P.A
                    bo = P.b

                    # osqp solver -> sd, p, λsd
                    Aex, bex = shift_to(Aes[i], bes[i], xt)
                    qq = [0; 0; 1.0]
                    AA = [Aex -bex; Ao -bo]
                    bb = [Aex * xt[1:2]; Ao * xt[1:2]]
                    # -AA * [p; sd] ≥ -bb
                    ret = solve_qp(UseOSQPSolver(); A=-sparse(AA), l=-bb, q=qq, polish=true, verbose=false)

                    p = ret.x[1:2]
                    sd = ret.x[3]
                    λsd = ret.y

                    #Main.@infiltrate
                    Jlag_rows, Jlag_cols, Jlag_vals, Jlag_buf = get_Jlag[i]
                    Jlag_vals(Jlag_buf, xt, Ao, bo, λsd, p, sd)
                    Jlag = sparse(Jlag_rows, Jlag_cols, Jlag_buf, 6, 8)

                    @inbounds JJ[xt_inds, xt_inds] .+= Jlag[1:n_x, 1:n_x]
                end

            end
        end
        nothing
    end

    n_z = length(z)
    n_nom = length(λ_nom)

    return (; fill_F!,
        J_both=(Jnom_rows, Jnom_cols, get_J_both!),
        l,
        u,
        T,
        n_z,
        n_xu,
        n_x,
        n_nom,
        n_obs,
        ego_polys,
        sides_per_poly,
        p1_max,
        p2_min)
end

function visualize_differentiablecol(x0, T, ego_polys, obs_polys; fig=Figure(), ax=Axis(fig[1, 1], aspect=DataAspect()), θ=[], is_displaying=true)
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

function solve_differentiablecol(prob, x0, obs_polys; θ0=nothing, is_displaying=true)
    (; fill_F!, J_both, ego_polys, l, u, T, n_z, n_xu, n_x, n_nom, n_obs, sides_per_poly, p1_max, p2_min) = prob

    @assert length(obs_polys) == n_obs
    J_rows, J_cols, J_vals! = J_both
    nnz_total = length(J_rows)

    n = length(l)
    @assert n == n_z + n_nom

    if is_displaying
        (fig, update_fig) = visualize_quick(x0, T, ego_polys, obs_polys)
    end

    if isnothing(θ0)
        θ0 = zeros(n)
        #derivs_per_fv = Int(n_α / length(ego_polys))
        #derivs_per_sd = Int(n_β / (T * length(ego_polys) * n_obs))
        for t in 1:T
            θ0[(t-1)*n_xu+1:(t-1)*n_xu+n_x] = x0
        end
        #θ0[n_z+1:n_z+n_α] .= 1.0 / derivs_per_fv
        #θ0[n_z+n_α+1:n_z+n_α+n_β] .= 1.0 / derivs_per_sd
    end

    #JJ = J_example
    #J_col = JJ.colptr[1:end-1]
    #J_len = diff(JJ.colptr)
    #J_row = JJ.rowval
    #nnz_total = length(JJ.nzval)

    J_shape = sparse(J_rows, J_cols, Vector{Cdouble}(undef, nnz_total), n, n)
    J_col = J_shape.colptr[1:end-1]
    J_len = diff(J_shape.colptr)
    J_row = J_shape.rowval


    function F(n, θ, result)
        result .= 0.0
        @inbounds z = θ[1:n_z]
        @inbounds λ_nom = θ[n_z+1:n_z+n_nom]
        fill_F!(result, z, x0, obs_polys, λ_nom)

        if is_displaying
            update_fig(θ)
        end

        Cint(0)
    end
    function J(n, nnz, θ, col, len, row, data)
        @assert nnz == nnz_total
        data .= 0.0
        @inbounds z = θ[1:n_z]
        @inbounds λ_nom = θ[n_z+1:n_z+n_nom]
        J_vals!(data, z, x0, obs_polys, λ_nom)
        col .= J_col
        len .= J_len
        row .= J_row
        #data .= JJ.nzval
        Cint(0)
    end


    buf = zeros(n)
    buf2 = zeros(n)
    Jbuf = zeros(nnz_total)

    w = randn(length(θ0))
    w = copy(θ0)

    F(n, w, buf)
    J(n, nnz_total, w, zero(J_col), zero(J_len), zero(J_row), Jbuf)

    #Jrows, Jcols, _ = findnz(J_example)
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

    if is_displaying
        display(fig)
    end

    @inbounds z = @view(θ[1:n_z])

    (; status, info, θ, z, fres)
end




