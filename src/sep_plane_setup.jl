function shift_to(V, x::AbstractArray{T}) where {T}
    p = x[1:2]
    θ = x[3]
    R = [cos(θ) sin(θ)
        -sin(θ) cos(θ)]
    Vx = map(V) do v
        R * v .+ p
    end
    Vx
end

# given a hyperplane a1*x + a2*y + b = 0, return line segment inside the limits
function get_line_segment(a1, a2, b; lim_x=[-10,10], lim_y=[-10,10])
    if a1 == 0
        xs = lim_x
        ys = [-b/a2, -b/a2]
    else
        ys = lim_y
        xs = (-a2 * ys .-b) ./ a1
    end
    return xs, ys
end

# deprecated
function g_col_sps(z, T, Vos, Ve; n_xu=9, n_sps=12)
    cons_sps = Num[]
    l_sps = Float64[]
    u_sps = Float64[]
    for t in 1:T
        xt = @view(z[(t-1)*n_xu+1:(t-1)*n_xu+3])
        Vex = shift_to(Ve, xt)
        m = length(Vex)
        for (e, V) in enumerate(Vos)
            ate = @view(z[n_xu*T+(t-1)*n_sps+(e-1)*3+1:n_xu*T+(t-1)*n_sps+(e-1)*3+2])
            bte = z[n_xu*T+(t-1)*n_sps+(e-1)*3+3]
            for i in 1:m # assuming number of vertices are equal between ego and obstacles
                push!(cons_sps, ate' * Vex[i] + bte)
                push!(l_sps, 0.0)
                push!(u_sps, Inf)
            end
            for i in 1:m
                push!(cons_sps, -ate' * V[i, :] - bte)
                push!(l_sps, 0.0)
                push!(u_sps, Inf)
            end
            push!(cons_sps, ate' * ate - 0.5)
            push!(l_sps, 0.0)
            push!(u_sps, Inf)
        end
    end
    (cons_sps, l_sps, u_sps)
end

function get_sps_cons(z, T, ego_polys, obs_polys, n_xu, n_per_col)
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
        xt = @view(z[(t-1)*n_xu+1:(t-1)*n_xu+3])
        yt = @view(z[T*n_xu+(t-1)*n_per_t+1:T*n_xu+(t-1)*n_per_t+n_per_t])

        for (i, Pi) in enumerate(ego_polys)
            yti = @view(yt[(i-1)*n_per_ego+1:(i-1)*n_per_ego+n_per_ego])
            Ve = Pi.V
            Vex = shift_to(Ve, xt)

            for (k, Pk) in enumerate(obs_polys)
                # ate = @view(yti[(k-1)*n_per_col+1:(k-1)*n_per_col+2])
                ate = yti[(k-1)*n_per_col+1]
                bte = yti[(k-1)*n_per_col+2]
                norm_unit = [cos(ate), sin(ate)]
                Vo = Pk.V

                for j in 1:sides_per_ego
                    push!(cons_sps, norm_unit' * Vex[j] + bte) # >= 0 
                    push!(l_sps, 0.0)
                    push!(u_sps, Inf)
                end

                for j in 1:sides_per_obs
                    push!(cons_sps, -norm_unit' * Vo[j] - bte) # >= 0 
                    push!(l_sps, 0.0)
                    push!(u_sps, Inf)
                end
            end
        end
    end
    (cons_sps, l_sps, u_sps)
end



function setup_sep_planes(
    ego_polys,
    obs_polys;
    T=1,
    dt=0.2,
    Rf=1e-3 * I(3),
    Qf=1e-3 * I(2),
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
    n_per_col = 2 # for every collision, we need two parameters, α and b, to represent a hyperplane cos(α)*x + sin(α)*y + b = 0 
    n_per_ego = n_per_col * n_obs
    n_per_t = n_per_col * n_obs * n_ego

    z = Symbolics.@variables(z[1:n_xu*T+n_per_t*T])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:xdim])[1] |> Symbolics.scalarize

    cost_nom = f(z, T, Rf, Qf)
    cons_dyn = g_dyn(z, x0, T, dt)
    cons_env = g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)

    #Ve = ego_polys[1].V
    #Vos = map(obs_polys) do P
    #    hcat(P.V...)' |> collect
    #end

    #cons_sps2, l_sps2, u_sps2 = g_col_sps(z, T, Vos, Ve)
    cons_sps, l_sps, u_sps = get_sps_cons(z, T, ego_polys, obs_polys, n_xu, n_per_col)

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
        p2_min,
        n_xu,
        obs_polys,
        n_per_col,
        n_per_ego,
        n_per_t
    )
end

function visualize_sep_planes(x0, T, ego_polys, obs_polys; n_per_col=3, fig=Figure(), ax=Axis(fig[1, 1], aspect=DataAspect()), θ=[], is_displaying=true, is_showing_sep_plane=true)
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    n_per_t = n_per_col * n_obs * n_ego
    n_xu = 9

    Vos = map(obs_polys) do P
        hcat(P.V...)' |> collect
    end

    xxts = Dict()
    abts = Dict()
    
    # limit of x and y to draw the hyperplane
    lim_x = [-10, 10]
    lim_y = [-10, 10]
    
    for i in 1:length(ego_polys)
        xx = x0[1:3]
        Aeb = shift_to(ego_polys[i].A, ego_polys[i].b, xx)
        self_poly = ConvexPolygon2D(Aeb[1], Aeb[2])
        Vex = self_poly.V

        plot!(ax, self_poly; color=:blue)
        for t in 1:T-1# 5:5:T-1
            xxts[i, t] = Observable(x0[1:3])
            Aeb = @lift(shift_to(ego_polys[i].A, ego_polys[i].b, $(xxts[i, t])))
            self_poly = @lift(ConvexPolygon2D($(Aeb)[1], $(Aeb)[2]))
            plot!(ax, self_poly; color=:blue, linestyle=:dash)

            if is_showing_sep_plane
                for (e, V) in enumerate(Vos)
                    abts[i, t, e] = Observable([0, -1.0])
                    # a1*x + a2*y + b = 0
                    a1 = @lift(cos($(abts[i, t, e])[1]))
                    a2 = @lift(sin($(abts[i, t, e])[1]))
                    b = @lift($(abts[i, t, e])[2])
                    line_seg = @lift(get_line_segment($a1, $a2, $b; lim_x=lim_x, lim_y=lim_y))
                    xs = @lift(($line_seg)[1])
                    ys = @lift(($line_seg)[2])
                    # GLMakie.lines!(ax, xs, ys; color=:green, linewidth=3)
                    # limits!(ax, lim_x, lim_y)
                end
            end
        end
        t = T
        xxts[i, t] = Observable(x0[1:3])
        Aeb = @lift(shift_to(ego_polys[i].A, ego_polys[i].b, $(xxts[i, t])))
        self_poly = @lift(ConvexPolygon2D($(Aeb)[1], $(Aeb)[2]))
        plot!(ax, self_poly; color=:blue, linewidth=3)

        if is_showing_sep_plane
            for e in 1:n_obs
                abts[i, t, e] = Observable([0, -1.0])
                # a1*x + a2*y + b = 0
                a1 = @lift(cos($(abts[i, t, e])[1]))
                a2 = @lift(sin($(abts[i, t, e])[1]))
                b = @lift($(abts[i, t, e])[2])
                line_seg = @lift(get_line_segment($a1, $a2, $b; lim_x=lim_x, lim_y=lim_y))
                xs = @lift(($line_seg)[1])
                ys = @lift(($line_seg)[2])
                GLMakie.lines!(ax, xs, ys; color=:green, linewidth=3)
                limits!(ax, lim_x, lim_y)
            end
        end
    end

    colors = [:red for _ in 1:n_obs]
    for (P, c) in zip(obs_polys, colors)
        plot!(ax, P; color=c)
    end

    function update_fig(θ)
        for i in 1:n_ego
            for t in 1:T
                xxts[i, t][] = copy(θ[(t-1)*n_xu+1:(t-1)*n_xu+6])
                if is_showing_sep_plane
                    for k in 1:n_obs
                        abts[i, t, k][] = copy(θ[n_xu*T+(t-1)*n_per_t+(k-1)*n_per_col+1:n_xu*T+(t-1)*n_per_t+(k-1)*n_per_col+2])
                    end
                end
            end
            #@infiltrate
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

function solve_prob_sep_planes(prob, x0; θ0=nothing, is_displaying=false, sleep_duration=0.0, is_recording=false)
    (; F_both!, J_both, l, u, T, n_z, n_nom, ego_polys, p1_max, p2_min, n_xu, obs_polys, n_per_col) = prob


    J_rows, J_cols, J_vals! = J_both
    nnz_total = length(J_rows)
    n = length(l)
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    n_per_ego = n_per_col * n_obs
    n_per_t = n_per_col * n_obs * n_ego

    @assert n == n_z + n_nom "did you forget to update l/u"

    if is_displaying
        (fig, update_fig) = visualize_sep_planes(x0, T, ego_polys, obs_polys; n_per_col)
    end

    if isnothing(θ0)
        θ0 = zeros(n)
        for t in 1:T
            θ0[(t-1)*n_xu+1:(t-1)*n_xu+6] = x0
            yt = @view(θ0[T*n_xu+(t-1)*n_per_t+1:T*n_xu+(t-1)*n_per_t+n_per_t])
            for (i, Pi) in enumerate(ego_polys)
                yti = @view(yt[(i-1)*n_per_ego+1:(i-1)*n_per_ego+n_per_ego])
                for (k, Pk) in enumerate(obs_polys)
                    #Vo = hcat(Pk.V...)' |> collect
                    #Vo_center = sum(Vo; dims=1) ./ size(Vo, 1)
                    #a = x0[1:2] - Vo_center[1:2]
                    #z = (x0[1:2] + Vo_center[1:2]) / 2
                    #b = -a'z

                    a = 0
                    b = -1
                    yti[(k-1)*n_per_col+1] = a
                    yti[(k-1)*n_per_col+2] = b
                end
            end
        end
        #for t in 1:T
        #    for i in 1:n_ego
        #        xxts[i, t][] = copy(θ0[(t-1)*n_xu+1:(t-1)*n_xu+6])
        #        for k in 1:n_obs
        #            abts[i, t, k][] = copy(θ0[n_xu*T+(t-1)*n_per_t+(k-1)*3+1:n_xu*T+(t-1)*n_per_t+(k-1)*3+3])
        #        end
        #    end
        #end
    end

    J_shape = sparse(J_rows, J_cols, Vector{Cdouble}(undef, nnz_total), n, n)
    J_col = J_shape.colptr[1:end-1]
    J_len = diff(J_shape.colptr)
    J_row = J_shape.rowval

    θ_history = []

    function F(n, θ, result)
        result .= 0.0
        @inbounds z = θ[1:n_z]
        @inbounds λ_nom = θ[n_z+1:n_z+n_nom]
        F_both!(result, z, x0, λ_nom)

        if is_recording
            push!(θ_history, copy(θ))
        end

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

    (; status, info, θ, z, fres, θ_history)
end
