function verts_from(x::AbstractArray{T}, angles, lengths) where {T}
    m = length(angles)
    V = zeros(T, m, 2)

    A, b = poly_from(x, angles, lengths)
    for i in 1:m-1
        V[i, :] = -A[i:i+1, :] \ b[i:i+1]
    end
    V[m, :] = -A[[m, 1], :] \ b[[m, 1]]
    V
end

function verts_from(A::AbstractArray{T}, b) where {T}
    m = length(b)
    V = zeros(T, m, 2)

    for i in 1:m-1
        V[i, :] = -A[i:i+1, :] \ b[i:i+1]
    end
    V[m, :] = -A[[m, 1], :] \ b[[m, 1]]
    V
end

function g_col_sps(z, T, V1, V2, V3, Ae, be)
    cons = Num[]
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+3])
        Aex, bex = shift_to(Ae, be, xt)
        verts = verts_from(Aex, bex)
        #@infiltrate
        m = size(verts, 1)
        for (e, V) in enumerate([V1, V2, V3])
            ate = @view(z[6*T+(t-1)*9+(e-1)*3+1:6*T+(t-1)*9+(e-1)*3+2])
            bte = z[6*T+(t-1)*9+(e-1)*3+3]
            for i in 1:m
                push!(cons, ate' * verts[i, :] + bte)
                push!(cons, -ate' * V[i, :] - bte)
            end
            #push!(cons, 1.0-ate'*ate)
            push!(cons, ate' * ate - 0.5)
        end
    end
    cons
end

function setup_sep_planes(ego_polys;
    T=40,
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
    rng=MersenneTwister(420),
    sides_per_poly=4,
    derivs_per_sd=4,
    derivs_per_fv=4,
    N_polys=3,
    angles=0:2*π/sides_per_poly:(2π-0.01),
    lengths=0.5 * rand(rng, sides_per_poly) .+ 0.25)

    xdim = 6

    n_sps = length(N_polys) * 3 * T
    z = Symbolics.@variables(z[1:xdim*T+n_sps])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:xdim])[1] |> Symbolics.scalarize

    cost = f(z, T, R)
    cons_dyn = g_dyn(z, x0, T, dt, L)
    cons_env = g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)
    #@infiltrate

    cons_nom = [cons_dyn; cons_env]
    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize

    θ = [z; λ_nom]

    lag = cost - cons_nom' * λ_nom
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; cons_nom]
    F_nom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false), parallel=Symbolics.SerialForm())[2]

    l = [fill(-Inf, length(grad_lag)); fill(-Inf, length(cons_dyn)); zeros(length(cons_env)); zeros(n_sps)]
    u = [fill(+Inf, length(grad_lag)); fill(Inf, length(cons_dyn)); fill(Inf, length(cons_env)); fill(Inf, n_sps)]
    n = length(l)

    J_nom = Symbolics.sparsejacobian(F_nom, θ)
    (J_rows_nom, J_cols_nom, J_vals) = findnz(J_nom)
    J_vals_nom! = Symbolics.build_function(J_vals, z, x0, λ_nom; expression=Val(false), parallel=Symbolics.SerialForm())[2]
	function F_both!(F, z, x0, polys, λ_nom, λ_sps)
        F .= 0.0
        F_nom!(F, z, x0, λ_nom)

        for i in 1:N_ego_polys
            Ae = collect(ego_polys[i].A)
            be = ego_polys[i].b

            for t in 1:T
                # ???
                xt_inds = (t-1)*n_xu+1:(t-1)*n_xu+xdim
                @inbounds xt = z[xt_inds]
                Aex, bex = shift_to(Ae, be, xt)
                verts = verts_from(Aex, bex)

                for (e, V) in enumerate(polys)
                    a_ind = n_xu*T+(i-1)*n_sps*T+(t-1)*n_sps+(e-1)*3+1:n_xu*T+(i-1)*n_sps*T+(t-1)*n_sps+(e-1)*3+2
                    b_ind = n_xu*T+(i-1)*n_sps*T+(t-1)*n_sps+(e-1)*3+3
                    λ_ind = (i-1)*n_sps*T+(t-1)*n_sps+(e-1)*3+1:(i-1)*n_sps*T+(t-1)*n_sps+(e-1)*3+3

                    ate = @view(z[a_ind])
                    bte = @view(z[b_ind])

                    for i in 1:sides_per_poly

                        F[n_xu*T+(t-1)*n_sps+(e-1)*3+1] .+= λ_ind[i]
                        F[n_xu*T+(t-1)*n_sps+(e-1)*3+2] .+=
                        F[n_xu*T+(t-1)*n_sps+(e-1)*3+1] .+=

                        push!(cons, ate' * verts[i, :] + bte)
                        push!(cons, -ate' * V[i, :] - bte)
                    end
                    #push!(cons, 1.0-ate'*ate)
                    push!(cons, ate' * ate)
                end
            end
        end
        nothing
    end


    function J_both_vals!(J_vals, z_local, x0_local, V1, V2, V3, λ_nom_local)
        J_vals .= 0.0
        J_vals_nom!(J_vals, z_local, V1, V2, V3, x0_local, λ_nom_local)
        nothing
    end


    return (; F_both!,
        J_both=(J_rows_nom, J_cols_nom, J_both_vals!),
        l,
        u,
        T,
        n_z=length(z),
        n_nom=length(λ_nom),
        N_polys,
        ego_polys,
        p1_max,
        p2_min)
end


function solve_prob_sep_planes(prob, x0, P1, P2, P3; θ0=nothing)
    (; F_both!, J_both, l, u, T, n_z, n_nom, N_polys, ego_polys, p1_max, p2_min) = prob

    J_rows, J_cols, J_vals! = J_both
    nnz_total = length(J_rows)
    n = length(l)
    @assert n == n_z + n_nom

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())

    polys = [P1, P2, P3]
    @assert length(polys) == N_polys

    xxts = Dict()
    for i in 1:length(ego_polys)
        xx = x0[1:3]
        Aeb = shift_to(ego_polys[i].A, ego_polys[i].b, xx)
        self_poly = ConvexPolygon2D(Aeb[1], Aeb[2])
        plot!(ax, self_poly; color=:blue)
        for t in 5:5:T-1
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

    display(fig)

    if isnothing(θ0)
        θ0 = zeros(n)
        T = Int(n_z / (6 + 3 * 3))
        for t in 1:T
            θ0[(t-1)*6+1:(t-1)*6+3] = x0[1:3]
            for e in 1:3
                θ0[6*T+(t-1)*9+(e-1)*3+1:6*T+(t-1)*9+(e-1)*3+2] = [0.0, 1]
                θ0[6*T+(t-1)*9+(e-1)*3+3] = x0[2] - 2.0
            end
        end
    end

    V1 = hcat(P1.V...)' |> collect
    V2 = hcat(P2.V...)' |> collect
    V3 = hcat(P3.V...)' |> collect

    #@infiltrate

    J_shape = sparse(J_rows, J_cols, Vector{Cdouble}(undef, nnz_total), n, n)
    J_col = J_shape.colptr[1:end-1]
    J_len = diff(J_shape.colptr)
    J_row = J_shape.rowval

    function F(n, θ, result)
        result .= 0.0
        #@inbounds z = @view(θ[1:n_z])
        #@inbounds λ_nom = @view(θ[n_z+1:n_z+n_nom])
        #@inbounds λ_col = @view(θ[n_z+n_nom+1:n_z+n_nom+n_col])
        @inbounds z = θ[1:n_z]
        @inbounds λ_nom = θ[n_z+1:n_z+n_nom]
        F_both!(result, z, x0, V1, V2, V3, λ_nom)
        Cint(0)
    end
    function J(n, nnz, θ, col, len, row, data)
        @assert nnz == nnz_total
        data .= 0.0
        @inbounds z = θ[1:n_z]
        @inbounds λ_nom = θ[n_z+1:n_z+n_nom]
        J_vals!(data, z, x0, V1, V2, V3, λ_nom)
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

    #@infiltrate

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
        cumulative_iteration_limit=50_000,
        #convergence_tolerance=1e-8
        convergence_tolerance=5e-4
    )
    #@infiltrate status != PATHSolver.MCP_Solved
    @inbounds z = @view(θ[1:n_z])
    @inbounds λ_nom = @view(θ[n_z+1:n_z+n_nom])

    (; status, info, θ, z, λ_nom)
end
