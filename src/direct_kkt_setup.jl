function get_direct_kkt_cons(z, T, ego_polys, obs_polys, n_xu, n_per_col)
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    n_per_ego = n_per_col * n_obs
    n_per_t = n_per_col * n_obs * n_ego

    cons_kkt = Num[]
    l_kkt = Float64[]
    u_kkt = Float64[]

    for t in 1:T
        xt = @view(z[(t-1)*n_xu+1:(t-1)*n_xu+3])
        yt = @view(z[T*n_xu+(t-1)*n_per_t+1:T*n_xu+(t-1)*n_per_t+n_per_t])

        for (i, Pi) in enumerate(ego_polys)
            yti = @view(yt[(i-1)*n_per_ego+1:(i-1)*n_per_ego+n_per_ego])
            Ae = Pi.A
            be = Pi.b
            Aex, bex = shift_to(Ae, be, xt)

            for (k, Pk) in enumerate(obs_polys)
                xxt = @view(yti[(k-1)*n_per_col+1:(k-1)*n_per_col+3])
                λt = @view(yti[(k-1)*n_per_col+4:(k-1)*n_per_col+n_per_col])
                Ao = Matrix(Pk.A)
                bo = Pk.b

                # min x' q, s.t. A x >= b
                AA, bb, qq = gen_LP_data(Aex, bex, Ao, bo)

                push!(cons_kkt, λt' * (AA * xxt + bb)) # = 0 (1)
                push!(l_kkt, -Inf)
                push!(u_kkt, Inf)

                append!(cons_kkt, qq - AA' * λt) # = 0 (3)
                append!(l_kkt, fill(-Inf, 3))
                append!(u_kkt, fill(+Inf, 3))

                append!(cons_kkt, AA * xxt + bb) # >= 0  (m1 + m2)
                append!(l_kkt, zeros(length(bb)))
                append!(u_kkt, fill(Inf, length(bb)))

                append!(cons_kkt, λt)  # >= 0 (m1 + m2)
                append!(l_kkt, zeros(length(bb)))
                append!(u_kkt, fill(Inf, length(bb)))

                push!(cons_kkt, xxt[3])
                append!(l_kkt, 0.0)
                append!(u_kkt, Inf)
            end
        end
    end
    (cons_kkt, l_kkt, u_kkt)
end

function setup_direct_kkt(
    ego_polys,
    obs_polys;
    T=1,
    dt=0.2,
    Rf=0.01 * I(3),
    Qf=0.01 * I(3),
    p1_max=500.0,
    p2_min=-500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=π / 4
)
    xdim = 6
    udim = 3
    n_xu = xdim + udim
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    sides_per_obs = length(obs_polys[1].b)
    sides_per_ego = length(ego_polys[1].b)
    n_per_col = sides_per_obs + sides_per_ego + 3
    n_per_ego = n_per_col * n_obs
    n_per_t = n_per_col * n_obs * n_ego

    z = Symbolics.@variables(z[1:n_xu*T+n_per_t*T])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:xdim])[1] |> Symbolics.scalarize


    cost_nom = f(z, T, Rf, Qf)
    cons_dyn = g_dyn(z, x0, T, dt)
    cons_env = g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)

    cons_kkt, l_kkt, u_kkt = get_direct_kkt_cons(z, T, ego_polys, obs_polys, n_xu, n_per_col)

    cons_nom = [cons_dyn; cons_env; cons_kkt]

    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize
    lag = cost_nom - cons_nom' * λ_nom
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; cons_nom]
    F_nom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false), parallel=Symbolics.SerialForm())[2]

    l = [fill(-Inf, length(grad_lag)); fill(-Inf, length(cons_dyn)); zeros(length(cons_env)); l_kkt]
    u = [fill(+Inf, length(grad_lag)); fill(Inf, length(cons_dyn)); fill(Inf, length(cons_env)); u_kkt]

    θ = [z; λ_nom]

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
        p1_max,
        p2_min,
        n_xu,
        n_per_col,
        obs_polys,
        ego_polys
    )
end


function solve_prob_direct_kkt(prob, x0; θ0=nothing)
    (; F_both!, J_both, l, u, T, n_z, n_nom, ego_polys, p1_max, p2_min, n_xu, n_per_col, obs_polys) = prob
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)

    J_rows, J_cols, J_vals! = J_both
    nnz_total = length(J_rows)
    n = length(l)

    @assert n == n_z + n_nom "did you forget to update l/u"

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())

    @assert length(obs_polys) == n_obs

    Vos = map(obs_polys) do P
        hcat(P.V...)' |> collect
    end

    xxts = Dict()
    abts = Dict()
    for i in 1:length(ego_polys)
        xx = x0[1:3]
        Aeb = shift_to(ego_polys[i].A, ego_polys[i].b, xx)
        self_poly = ConvexPolygon2D(Aeb[1], Aeb[2])

        plot!(ax, self_poly; color=:blue)
        for t in 1:T# 5:5:T-1
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

    display(fig)

    if isnothing(θ0)
        θ0 = zeros(n)
        for t in 1:T
            θ0[(t-1)*n_xu+1:(t-1)*n_xu+6] = x0
        end

        for i in 1:length(ego_polys)
            for t in 1:T
                xxts[i, t][] = copy(θ0[(t-1)*n_xu+1:(t-1)*n_xu+6])
            end
        end
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

        for i in 1:length(ego_polys)
            for t in 1:T
                xxts[i, t][] = copy(θ[(t-1)*n_xu+1:(t-1)*n_xu+6])
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

    buf = zeros(n)
    buf2 = zeros(n)
    Jbuf = zeros(nnz_total)

    F(n, θ0, buf)
    J(n, nnz_total, θ0, zero(J_col), zero(J_len), zero(J_row), Jbuf)



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
        #convergence_tolerance=1e-8
        convergence_tolerance=5e-4
    )


    @inbounds z = @view(θ[1:n_z])
    @inbounds λ_nom = @view(θ[n_z+1:n_z+n_nom])

    (; status, info, θ, z, λ_nom)
end
