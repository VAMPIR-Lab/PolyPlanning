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
    xprev = x0
    for t in 1:T
        xt = @view(x[(t-1)*6+1:(t-1)*6+4])
        ut = @view(x[(t-1)*6+5:(t-1)*6+6])
        append!(g, xt - kinematic_bicycle_dyn(x_prev, ut, dt))
        x_prev = xt
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
                    sds[t, e, ee] = sd
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

function setup()
    T = 3
    z = Symbolics.@variables(z[1:6*T])[1] |> Symbolics.scalarize
    angles = [0.0, π/2, π, 3*π/2]
    lengths = [0.5, 1.0, 0.5, 0.5]
    P1 = ConvexPolygon2D([[0.0, 0.0], [1,1], [1,2], [0,3]])
    sds = g_col_all(z, T, [P1], angles, lengths)
    (sds, z)
end
