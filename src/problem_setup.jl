function f(z, T, goal_dir, R)
    cost = 0.0
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+4])
        ut = @view(z[(t-1)*6+5:(t-1)*6+6])
        cost += ut'*R*ut - goal_dir'*xt[1:2]
    end
    cost
end

function pointmass_dyn(x, u, dt)
    p1, p2, v1, v2 = x
    a1, a2 = u
    x + dt * [v1 + dt/2 * a1, v2 + dt/2 * a2, a1, a2]
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
        append!(g, xt - kinematic_bicycle_dyn(x_prev, ut, dt, L))
        x_prev = xt
    end
    g
end

function g_env(z, T, p1_max)
    g = Num[]
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+4])
        append!(g, [xt[1] - p1_max, p1_max-xt[1], xt[2]])
    end
    g
end


function poly_from(x, angles, lengths)
    m = length(angles)
    A = zeros(Num, m, 2)
    b = zeros(Num, m)
    p1, p2, θ, v = x
    p = [p1,p2]
    for i in 1:m
        θi = θ + angles[i]
        ai = [cos(θi), sin(θi)]
        bi = lengths[i] - ai'*p
        A[i,:] += ai
        b[i] += bi
    end
    A, b
end

function gen_lmcp_data(A1,b1,A2,b2) 
    m1 = length(b1)
    m2 = length(b2)
    M = [zeros(Num, m1,m1) A1 zeros(Num, m1,m2) -ones(Num, m1,1);
         -A1'  zeros(2, 2) -A2' zeros(2,1);
         zeros(Num, m2, m1) A2 zeros(Num, m2,m2) zeros(Num, m2,1);
         ones(Num, 1,m1) zeros(Num, 1,2) zeros(Num, 1,m2) zeros(Num, 1,1)]
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
            M, q, l, u = gen_lmcp_data(A,b,poly.A,poly.b)
            m1 = length(b)
            m2 = length(poly.b)
            Itr = Iterators.product([ 1 ≤ i ≤ m1 || m1+3 ≤ i ≤ m1+m2+2 ? [true, false] : [true,] for i in 1:m1+m2+3]...)
            for (ee,assignment) in enumerate(Itr)
                assignment = collect(assignment)
                zz = zeros(Num, m1+m2+3)
                try
                    zz[assignment] = M[assignment, assignment] \ q[assignment]
                    yy = zz[1:m1] 
                    xx = zz[m1+1:m1+2]
                    sd = yy'*(A*xx+b)
                    sds[t, e, assignment] = sd
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
            M, q, l, u = gen_lmcp_data(A,b,poly.A,poly.b)
            m1 = length(b)
            m2 = length(poly.b)
            ret = solve_lmcp(UsePATHSolver(), M, q, l, u) 
            z = ret.z
            I1 = z .≥ (l .+ 1e-4) .&& (M*z+q .< 1e-4)
            I2 = z .< (l .+ 1e-4) .&& (M*z+q .< 1e-4)
            I3 = z .< (l .+ 1e-4) .&& (M*z+q .≥ 1e-4)

            y = ret.z[1:m1]
            x = ret.z[m1+1:m1+2]
            sd = y'*(A1*x+b1)
            
            if sum(I2) > 1
                error("More derivatives than accounted for")
            elseif sum(I2) == 1
                assignment1 = I1 .|| I2
                assignment2 = I1
            else
                assignment1 = I1
                assignment2 = I1
            end
            IDs[t,e] = (assignment1, assignment2)
            SDs[t,e] = sd
        end
    end
end

function setup()
    T = 1
    dt = 1.0
    L = 1.0
    goal_dir = [0, -1.0]
    R = 0.01*I(2,2)
    p1_max = 1.0
    rng = MersenneTwister(420)
    
    P1 = ConvexPolygon2D([randn(rng, 2) for _ in 1:5] .+ [-3.0,0])
    P2 = ConvexPolygon2D([randn(rng, 2) for _ in 1:5] .+ [3.0,0])
    P3 = ConvexPolygon2D([randn(rng, 2) for _ in 1:5] .+ [0.0,-1])
    avoid_polys = [P1,P2,P3]
    N_polys = length(avoid_polys)
    sds = g_col_all(z, T, avoid_polys, angles, lengths)

    z = Symbolics.@variables(z[1:6*T+4])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:4])[1] |> Symbolics.scalarize

    angles = [0.0, π/2, π, 3*π/2]
    lengths = [0.5, 1.0, 0.5, 0.5]
    sd_funcs = Dict()
    sd_grads = Dict()
    sd_hesss = Dict()
    for (k,sd) in sds  
        grad = Symbolics.gradient(sd, z)
        hess = Symbolics.sparsejacobian(grad, z)
        (H_rows, H_cols, H_vals) = findnz(hess)

        get_sd = Symbolics.build_function(sd, z; expression=Val(false))
        get_grad = Symbolics.build_function(grad, z; expression=Val(false))[1] # not in place version for this one
        get_H_vals = Symbolics.build_function(H_vals, z; expression=Val(false))[1] # not in place version for this one

        sd_funcs[k] = get_sd
        sd_grads[k] = get_grad
        sd_hesss[k] = (H_rows, H_cols, get_H_vals)
    end

    cost = f(z, T, goal_dir, R)
    cons_nom = [g_dyn(z, x0, T, dt, L); g_env(z, T, p1_max)]
    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize    

    lag = cost - cons_nom'*λ_nom
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; cons_nom; zeros(Num, T*N_polys)]
    J_nom = Symbolics.sparsejacobian(F_nom, z)
    (J_rows_nom, J_cols_nom, J_vals) = findnz(J_nom)

    F_nom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false))[2]
    J_vals_nom! = Symbolics.build_function(J_vals, z, x0, λ_nom; expression=Val(false))[2]

    function F_col!(F, z, λ_col)
        (sds, sd_ids) = get_sd_ids(z, T, avoid_polys, angles, lengths)
        for t in 1:T
            for e in 1:N_polys
                assignments = sd_ids[t,e]
                F[1:6*T] -= sd_grads[t,e,assignments[1]] * λ_col[(t-1)*2*N_polys+(e-1)*2+1]
                F[1:6*T] -= sd_grads[t,e,assignments[2]] * λ_col[(t-1)*2*N_polys+(e-1)*2+2]
                F[6*T+length(cons_nom)+(t-1)*N_polys+e] += sds[t,e]
            end
        end
    end
    function J_col!(F,z,λ_col)
        #todo
    end

end
