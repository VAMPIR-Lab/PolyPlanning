function f(z, T, goal_dir, R)
    cost = 0.0
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+4])
        ut = @view(z[(t-1)*6+5:(t-1)*6+6])
        #cost += ut'*R*ut - goal_dir'*xt[1:2]
        #cost += ut'*R*ut + xt[1:2]'*xt[1:2]
        cost += xt[1:2]'*xt[1:2]
    end
    cost
end

function pointmass_dyn(x, u, dt)
    p1, p2, v1, v2 = x
    a1, a2 = u
    x + [a1, a2, 0, 0]
    #x + dt * [v1 + dt/2 * a1, v2 + dt/2 * a2, a1, a2]
end

function kinematic_bicycle_dyn(x, u, dt, L)
    p1, p2, θ, v = x
    δ, a = u
    x + dt * [cos(θ)*v, sin(θ)*v, tan(δ)*v / L, a]
end

function g_dyn(z, x0, T, dt, L)
    g = Num[] 
    x_prev = x0
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+4])
        ut = @view(z[(t-1)*6+5:(t-1)*6+6])
        #append!(g, xt - kinematic_bicycle_dyn(x_prev, ut, dt, L))
        append!(g, xt - pointmass_dyn(x_prev, ut, dt))
        x_prev = xt
    end
    g
end

function g_env(z, T, p1_max)
    g = Num[]
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+4])
        append!(g, [p1_max + xt[1], p1_max-xt[1], xt[2]])
    end
    g
end


function poly_from(x::AbstractArray{T}, angles, lengths) where T
    m = length(angles)
    A = zeros(T, m, 2)
    b = zeros(T, m)
    p1, p2, θ, v = x
    p = [p1,p2]
    for i in 1:m
        #θi = θ + angles[i]
        θi = angles[i]
        ai = [cos(θi), sin(θi)]
        bi = lengths[i] - ai'*p
        A[i,:] += ai
        b[i] += bi
    end
    A, b
end

function gen_LP_data(A1::AbstractArray{T},b1,A2,b2) where T
    m1 = length(b1)
    m2 = length(b2)
    A = [[A1; A2] ones(T, m1+m2)]
    b = [b1; b2]
    q = [0, 0, 1.0]
    (A, b, q)
end

function get_LP_signed_distance(P1, P2)
    AA, bb, qq = gen_LP_data(P1.A,P1.b,P2.A,P2.b)
    m1 = length(P1.b)
    m2 = length(P2.b)

    ret = solve_qp(UseOSQPSolver(); A=sparse(AA), l=-bb, q=qq, polish=true, verbose=false)
    primals = ret.x
    duals = -ret.y
    cons = AA*primals+bb

    I1 = duals .≥ 1e-4 .&& cons .< 1e-4
    I2 = duals .< 1e-4 .&& cons .< 1e-4
    I3 = duals .< 1e-4 .&& cons .≥ 1e-4

    active_inds = collect(1:m1+m2)[I1]
    M = [zeros(3, 3) -AA';
         AA zeros(m1+m2,m1+m2)]
    r = [qq; bb]
    AA_active = collect(AA[active_inds,:])
    bb_active = collect(bb[active_inds])
    if length(active_inds) ≥ 3
        zz = -(AA_active'*AA_active)\AA_active'*bb_active
        dsd = ((AA_active'*AA_active)\AA_active')*AA_active
    else
        zz = -AA_active'*((AA_active*AA_active')\bb_active)
        dsd = AA_active'*((AA_active*AA_active')\AA_active)
    end
    sd2 = zz[3]
    dsd2 = dsd[3,:]

    sd = primals[3]
    (; sd, sd2, dsd2, I1, I2, I3, AA_active, bb_active, active_inds)
end

function gen_lmcp_data(A1::AbstractArray{T},b1,A2,b2) where T
    m1 = length(b1)
    m2 = length(b2)
    M = [zeros(T, m1,m1) A1 zeros(T, m1,m2) -ones(T, m1,1);
         -A1'  zeros(2, 2) -A2' zeros(2,1);
         zeros(T, m2, m1) A2 zeros(T, m2,m2) zeros(T, m2,1);
         ones(T, 1,m1) zeros(T, 1,2) zeros(T, 1,m2) zeros(T, 1,1)]
    q = [b1; zeros(2); b2; -1]
    l = [zeros(m1); fill(-Inf, 2); zeros(m2); fill(-Inf, 1)]
    u = fill(Inf, m1+m2+3)
    M, q, l, u
end

function g_col_all(z, T, avoid_polys, angles, lengths) 
    sds = Dict()
    num_polys = length(avoid_polys)
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+4])
        A,b = poly_from(xt, angles, lengths)
        for (e,poly) in enumerate(avoid_polys)
            AA, bb, qq = gen_LP_data(A,b,poly.A,poly.b)
            m1 = length(b)
            m2 = length(poly.b)
            M = [zeros(Num, 3, 3) -AA';
                 AA zeros(Num, m1+m2,m1+m2)]
            r = [qq; bb]
            
            all_active_inds = collect(1:m1+m2)
            Itr = powerset(all_active_inds) |> collect
            #Itr = Iterators.product([ 1 ≤ i ≤ m1 || m1+3 ≤ i ≤ m1+m2+2 ? [true, false] : [true,] for i in 1:m1+m2+3]...)
            for (ee,active_inds) in enumerate(Itr)
                #assignment = [[true, true, true];[i ∈ active_inds for i in 1:m1+m2]]
                assignment = [i ∈ active_inds for i in 1:m1+m2]
                #zz = zeros(Num, m1+m2+3)
                try
                    AA_active = collect(AA[active_inds,:])
                    bb_active = collect(bb[active_inds])
                    if length(active_inds) ≥ 3
                        zz = -(AA_active'*AA_active)\AA_active'*bb_active
                    else
                        zz = -AA_active'*((AA_active*AA_active')\bb_active)
                    end
                    #zz[assignment] = -M[assignment, assignment] \ r[assignment]
                    sd = zz[3]
                    #xx = zz[1:2]
                    sds[t, e, active_inds] = sd
                catch err
                    if err isa LinearAlgebra.SingularException
                        continue
                    else
                        @infiltrate
                    end
                end
            end
        end
    end
    sds
end

function get_sd_ids(z, T, avoid_polys, angles, lengths)
    num_polys = length(avoid_polys) 
    IDs = Dict()
    SDs = Dict()
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+4])
        A,b = poly_from(xt, angles, lengths)
        for (e,poly) in enumerate(avoid_polys)
            AA, bb, qq = gen_LP_data(A,b,poly.A,poly.b)
            m1 = length(b)
            m2 = length(poly.b)

            ret = solve_qp(UseOSQPSolver(); A=sparse(AA), l=-bb, q=qq, polish=true, verbose=false)
            primals = ret.x
            duals = -ret.y

            cons = AA*primals+bb
            I1 = duals .≥ 1e-4 .&& cons .< 1e-4
            I2 = duals .< 1e-4 .&& cons .< 1e-4
            I3 = duals .< 1e-4 .&& cons .≥ 1e-4
            sd = primals[3]
           
            all_inds = collect(1:m1+m2)
            if sum(I2) > 1
                error("More derivatives than accounted for! |I2| = $(sum(I2)). λ = $duals. cons = $cons")
            elseif sum(I2) == 1
                @info "two derivatives"
                assignment1 = all_inds[I1 .|| I2]
                assignment2 = all_inds[I1]
            else
                assignment1 = all_inds[I1]
                assignment2 = all_inds[I1]
            end
            IDs[t,e] = (assignment1, assignment2)
            SDs[t,e] = sd
        end
    end
    SDs, IDs
end

#function g_col_all(z, T, avoid_polys, angles, lengths) 
#    sds = Dict()
#    num_polys = length(avoid_polys)
#    for t in 1:T
#        xt = @view(z[(t-1)*6+1:(t-1)*6+4])
#        A,b = poly_from(xt, angles, lengths)
#        for (e,poly) in enumerate(avoid_polys)
#            M, q, l, u = gen_lmcp_data(A,b,poly.A,poly.b)
#            m1 = length(b)
#            m2 = length(poly.b)
#            Itr = Iterators.product([ 1 ≤ i ≤ m1 || m1+3 ≤ i ≤ m1+m2+2 ? [true, false] : [true,] for i in 1:m1+m2+3]...)
#            for (ee,assignment) in enumerate(Itr)
#                assignment = collect(assignment)
#                zz = zeros(Num, m1+m2+3)
#                try
#                    zz[assignment] = M[assignment, assignment] \ q[assignment]
#                    yy = zz[1:m1] 
#                    xx = zz[m1+1:m1+2]
#                    sd = yy'*(A*xx+b)
#                    sds[t, e, assignment] = sd
#                catch err
#                    if err isa LinearAlgebra.SingularException
#                        continue
#                    else
#                        @infiltrate
#                    end
#                end
#            end
#        end
#    end
#    sds
#end
#
#function get_sd_ids(z, T, avoid_polys, angles, lengths)
#    num_polys = length(avoid_polys) 
#    IDs = Dict()
#    SDs = Dict()
#    for t in 1:T
#        xt = @view(z[(t-1)*6+1:(t-1)*6+4])
#        A,b = poly_from(xt, angles, lengths)
#        for (e,poly) in enumerate(avoid_polys)
#            M, q, l, u = gen_lmcp_data(A,b,poly.A,poly.b)
#            m1 = length(b)
#            m2 = length(poly.b)
#            ret = solve_lmcp(UsePATHSolver(), M, q, l, u) 
#            z = ret.z
#            I1 = z .≥ (l .+ 1e-4) .&& (M*z+q .< 1e-4)
#            I2 = z .< (l .+ 1e-4) .&& (M*z+q .< 1e-4)
#            I3 = z .< (l .+ 1e-4) .&& (M*z+q .≥ 1e-4)
#
#            y = ret.z[1:m1]
#            x = ret.z[m1+1:m1+2]
#            sd = y'*(A*x+b)
#            
#            if sum(I2) > 1
#                error("More derivatives than accounted for")
#            elseif sum(I2) == 1
#                assignment1 = I1 .|| I2
#                assignment2 = I1
#            else
#                assignment1 = I1
#                assignment2 = I1
#            end
#            IDs[t,e] = (assignment1, assignment2)
#            SDs[t,e] = sd
#        end
#    end
#    SDs, IDs
#end

function setup(;
                T = 1,
                dt = 1.0,
                L = 1.0,
                goal_dir = [0, -10.0],
                R = 0.01*I(2),
                p1_max = 5.0,
                rng = MersenneTwister(421),
                angles = 0:2*π/5:(2π-0.01),
                lengths = 0.1*(1:5).+1)
    
    P1 = ConvexPolygon2D([randn(rng, 2) + [-3,0] for _ in 1:8])
    P1 = ConvexPolygon2D([[-2,-5],[-4,2.0], [-6,-5]])
    P2 = ConvexPolygon2D([randn(rng, 2) + [ 3,0] for _ in 1:8])
    P3 = ConvexPolygon2D([randn(rng, 2) + [0,-1] for _ in 1:8])
    avoid_polys = [P1,]
    N_polys = length(avoid_polys)

    z = Symbolics.@variables(z[1:6*T])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:4])[1] |> Symbolics.scalarize
    sds = g_col_all(z, T, avoid_polys, angles, lengths)
    
    cost = f(z, T, goal_dir, R)
    cons_dyn = g_dyn(z, x0, T, dt, L)
    cons_env = g_env(z, T, p1_max)
    cons_nom = [cons_dyn; cons_env]
    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize    
    λ_col = Symbolics.@variables(λ_col[1:T*2*N_polys])[1] |> Symbolics.scalarize

    θ = [z; λ_nom; λ_col]

    lag = cost - cons_nom'*λ_nom
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; cons_nom; zeros(Num, T*2*N_polys)]
    
    l = [fill(-Inf, length(grad_lag)); fill(-Inf, length(cons_dyn)); zeros(length(cons_env)); zeros(T*2*N_polys)]
    u = [fill(+Inf, length(grad_lag)); fill(Inf, length(cons_dyn)); fill(Inf, length(cons_env)); fill(Inf, (T*2*N_polys))]
    n = length(l)

    J_nom = Symbolics.sparsejacobian(F_nom, θ)
    (J_rows_nom, J_cols_nom, J_vals) = findnz(J_nom)
    J_nom = sparse(J_rows_nom, J_cols_nom, J_vals, n, n)

    F_nom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false))[2]

    sd_F_funcs = Dict()
    sd_Js = Dict()
    sd_J_funcs = Dict()

    
    J_pattern = collect(zip(J_rows_nom, J_cols_nom))
    
    num_sds = length(sds)

    i = 0
    for (k,sd) in sds  
        t = k[1]
        e = k[2]
        λ1 = λ_col[(t-1)*2*N_polys + (e-1)*2+1]
        λ2 = λ_col[(t-1)*2*N_polys + (e-1)*2+2]

        lag1 = Symbolics.gradient(-sd*λ1, z) |> simplify
        lag2 = Symbolics.gradient(-sd*λ2, z) |> simplify
        cons1 = zeros(Num, T*2*N_polys)
        cons2 = zeros(Num, T*2*N_polys)
        cons1[(t-1)*2*N_polys+(e-1)*2+1] = sd
        cons2[(t-1)*2*N_polys+(e-1)*2+2] = sd

        F1 = [lag1; zeros(Num, length(cons_nom)); cons1]
        F2 = [lag2; zeros(Num, length(cons_nom)); cons2]
    
        J1 = Symbolics.sparsejacobian(F1, θ)
        J2 = Symbolics.sparsejacobian(F2, θ)
        (J1_rows, J1_cols, J1_vals) = findnz(J1)
        (J2_rows, J2_cols, J2_vals) = findnz(J2)
        J1 = sparse(J1_rows, J1_cols, J1_vals, n, n)
        J2 = sparse(J2_rows, J2_cols, J2_vals, n, n)
        sd_Js[k] = (; loc_1 = J1, loc_2 = J2)

        J1_pattern = collect(zip(J1_rows, J1_cols))
        J2_pattern = collect(zip(J2_rows, J2_cols))
        union!(J_pattern, J1_pattern)
        union!(J_pattern, J2_pattern)

        get_F1 = Symbolics.build_function(F1, z, λ_col; expression=Val(false))[2]
        get_F2 = Symbolics.build_function(F2, z, λ_col; expression=Val(false))[2]
        
        sd_F_funcs[k] = (; loc_1=get_F1, loc_2=get_F2)
        i+=1
        @info "$i / $num_sds"
    end
    J_both_rows, J_both_cols = map(x->getfield.(J_pattern, x), fieldnames(eltype(J_pattern)))
    zero_J_pattern = sparse(J_both_rows, J_both_cols, fill(2*eps(),length(J_both_cols)), n, n)
    i = 0
    for k in keys(sds)
        Js = sd_Js[k]
        J1 = Js.loc_1 + zero_J_pattern
        J2 = Js.loc_2 + zero_J_pattern
        (J1_rows, J1_cols, J1_vals) = findnz(J1)
        J1_vals_perm = zeros(Num, length(J1_vals))
        for (e,(row,col)) in enumerate(zip(J1_rows, J1_cols))
            j = findfirst([ind == (row,col) for ind in J_pattern])
            J1_vals_perm[j] = J1_vals[e]
        end

        (J2_rows, J2_cols, J2_vals) = findnz(J2)
        J2_vals_perm = zeros(Num, length(J2_vals))
        for (e,(row,col)) in enumerate(zip(J2_rows, J2_cols))
            j = findfirst([ind == (row,col) for ind in J_pattern])
            J2_vals_perm[j] = J2_vals[e]
        end
        get_J1_vals! = Symbolics.build_function(J1_vals_perm, z, λ_col; expression=Val(false))[2]
        get_J2_vals! = Symbolics.build_function(J2_vals_perm, z, λ_col; expression=Val(false))[2]
        sd_J_funcs[k] = (; loc_1=get_J1_vals!, loc_2=get_J2_vals!)
        i+=1
        @info "$i / $num_sds"
    end
   
    J_nom += zero_J_pattern
    (J_nom_rows, J_nom_cols, J_vals) = findnz(J_nom)
    J_vals_perm = zeros(Num, length(J_vals))
    for (e,(row,col)) in enumerate(zip(J_nom_rows, J_nom_cols))
        j = findfirst([ind == (row,col) for ind in J_pattern])
        J_vals_perm[j] = J_vals[e]
    end
    J_vals_nom! = Symbolics.build_function(J_vals_perm, z, x0, λ_nom; expression=Val(false))[2]

    F_nom_buf = zeros(length(F_nom))
    F_col1_buf = zeros(length(F_nom))
    F_col2_buf = zeros(length(F_nom))
     
    function F_both!(F, z_local, x0_local, λ_nom_local, λ_col_local)
        F .= 0.0
        F_nom_buf .= 0.0
        F_col1_buf .= 0.0
        F_col2_buf .= 0.0

        F_nom!(F_nom_buf, z_local, x0_local, λ_nom_local)
        F .+= F_nom_buf
        sds, sd_ids = get_sd_ids(z_local, T, avoid_polys, angles, lengths)
        for t in 1:T
            for e in 1:N_polys
                assignments = sd_ids[t,e]
                F_func_1! = sd_F_funcs[t,e,assignments[1]].loc_1
                F_func_2! = sd_F_funcs[t,e,assignments[2]].loc_2
                F_func_1!(F_col1_buf, z_local, λ_col_local)
                F_func_2!(F_col2_buf, z_local, λ_col_local)
                F .+= F_col1_buf
                F .+= F_col2_buf
            end
        end
        nothing
    end

    J_nom_buf = zeros(length(J_vals))
    J_col1_buf = zeros(length(J_vals))
    J_col2_buf = zeros(length(J_vals))

    function J_both_vals!(J_vals, z_local, x0_local, λ_nom_local, λ_col_local)
        J_vals .= 0.0
        J_nom_buf .= 0.0
        J_col1_buf .= 0.0
        J_col2_buf .= 0.0

        J_vals_nom!(J_nom_buf, z_local, x0_local, λ_nom_local)
        J_vals .+= J_nom_buf
        sds, sd_ids = get_sd_ids(z_local, T, avoid_polys, angles, lengths)
        for t in 1:T
            for e  in 1:N_polys
                assignments = sd_ids[t,e]
                J_func_1! = sd_J_funcs[t,e,assignments[1]].loc_1
                J_func_2! = sd_J_funcs[t,e,assignments[2]].loc_2
                J_func_1!(J_col1_buf, z_local, λ_col_local)
                J_func_2!(J_col2_buf, z_local, λ_col_local)
                J_vals .+= J_col1_buf
                J_vals .+= J_col2_buf
            end
        end
        nothing
    end
    
    return (; F_both!, 
            J_both=(J_both_rows, J_both_cols, J_both_vals!), 
            l, 
            u, 
            n_z=length(z), 
            n_nom=length(λ_nom), 
            n_col=length(λ_col), 
            avoid_polys, 
            angles, 
            lengths,
            p1_max)
end

function solve_prob(prob, x0; θ0 = nothing)
    (; F_both!, J_both, l, u, n_z, n_nom, n_col, avoid_polys, angles, lengths, p1_max) = prob

    J_rows, J_cols, J_vals! = J_both
    nnz_total = length(J_rows)
    n = length(l)
    @assert n == n_z + n_nom + n_col

    fig = Figure()
    ax = Axis(fig[1,1], aspect=DataAspect())

    A,b = poly_from(x0, angles, lengths)
    self_poly = ConvexPolygon2D(A,b)

    plot!(ax, self_poly; color=:blue)
    for P in avoid_polys
        plot!(ax, P; color=:red)
    end
    lines!(ax, [-p1_max, p1_max],[0, 0], color=:black)
    lines!(ax, [-p1_max, -p1_max], [0, 10], color=:black)
    lines!(ax, [p1_max, p1_max], [0, 10], color=:black)
    display(fig)

    if isnothing(θ0) 
        θ0 = zeros(n)
        T = Int(n_z / 6)
        for t in 1:T
            θ0[(t-1)*6+1:(t-1)*6+4] = x0
        end
    end


    J_shape = sparse(J_rows, J_cols, Vector{Cdouble}(undef, nnz_total), n, n)
    J_col = J_shape.colptr[1:end-1]
    J_len = diff(J_shape.colptr)
    J_row = J_shape.rowval

    function F(n, θ, result)
        result .= 0.0
        @inbounds z = @view(θ[1:n_z])
        @inbounds λ_nom = @view(θ[n_z+1:n_z+n_nom])
        @inbounds λ_col = @view(θ[n_z+n_nom+1:n_z+n_nom+n_col])
        F_both!(result, z, x0, λ_nom, λ_col)
        Cint(0)
    end
    function J(n, nnz, θ, col, len, row, data)
        @assert nnz == nnz_total
        data .= 0.0
        @inbounds z = @view(θ[1:n_z])
        @inbounds λ_nom = @view(θ[n_z+1:n_z+n_nom])
        @inbounds λ_col = @view(θ[n_z+n_nom+1:n_z+n_nom+n_col])
        J_vals!(data, z, x0, λ_nom, λ_col)
        col .= J_col
        len .= J_len
        row .= J_row
        Cint(0)
    end

    PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")
    status, θ, info = PATHSolver.solve_mcp(
         F,
         J,
         l,
         u,
         θ0;
         silent=false,
         nnz=nnz_total,
         jacobian_structure_constant = true,
         output_linear_model = "no",
         preprocess = 0,
         output_warnings = "no",
         jacobian_data_contiguous = true,
         cumulative_iteration_limit = 50_000,
         convergence_tolerance=1e-8,
     )
    

    
    buf = zeros(n)
    buf2 = zeros(n)
    F(n, θ, buf)
    
    Jbuf = zeros(nnz_total)
    J(n, nnz_total, θ, zero(J_col), zero(J_len), zero(J_row), Jbuf)
    JJ = sparse(J_rows, J_cols, Jbuf, n, n) |> collect
    JJnum = zeros(n,n)
    JJact = zeros(n,n)

    for i in 1:n
        dθ = randn(n)
        dθ = dθ / norm(dθ)
        F(n,θ+1e-5*dθ, buf2)
        JJnum[:,i] = (buf2-buf) / 1e-5
        JJact[:,i] = JJ*dθ
    end

    sds = buf[end-n_col+1:end]
    A,b = poly_from(θ[1:4], angles, lengths)
    self_poly = ConvexPolygon2D(A,b)

    plot!(ax, self_poly; color=:blue, linestyle=:dash)
    display(fig)

    @infiltrate

    #@infiltrate status != PATHSolver.MCP_Solved
    @inbounds z = @view(θ[1:n_z])
    @inbounds λ_nom = @view(θ[n_z+1:n_z+n_nom])
    @inbounds λ_col = @view(θ[n_z+n_nom+1:n_z+n_nom+n_col])
    
    (; status, info, θ, z, λ_nom, λ_col, sds, JJ)
end
