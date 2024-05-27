# for given poly represented by A and b, translate it by x[1:2], and rotate it clockwise around point c by x[3]
function shift_to_around_point(P, x; c = 0)
    A = P.A
    b = P.b
    V = P.V
    if c==0
        c = sum(V) / length(V)
    end
    At, bt = PolyPlanning.shift_to(A, b, [-c[1], -c[2], 0.])
    At, bt = PolyPlanning.shift_to(At, bt, [0., 0., x[3]])
    At, bt = PolyPlanning.shift_to(At, bt, [c[1], c[2], 0.])
    At, bt = PolyPlanning.shift_to(At, bt, [x[1], x[2], 0.])
    ConvexPolygon2D(At, bt)
end

# for given poly P and two state x1, x2, generate polys interpolated between these two states
function interpolation_polys(P, x1, x2; c = 0, num = 10)
    if c==0
        c = sum(P.V) / length(P.V)
    end
    polys = []
    for t in 0:1/(num-1):1
        x = x1 * (1-t)+ x2 * t
        poly = PolyPlanning.shift_to_around_point(P, x; c=c)
        append!(polys, [poly])
    end
    polys
end

function get_LP_data_point_shrunk(P1, P2; c1=0, c2=0)
    if c1 == 0
        c1 = sum(P1.V) / length(P1.V)
    end
    if c2 == 0
        c2 = sum(P2.V) / length(P2.V)
    end

    dim = length(c1)
    P = spzeros(dim+1, dim+1)
    q = [zeros(dim); 1]
    A1 = P1.A
    b1 = P1.b
    A2 = P2.A
    b2 = P2.b
    A = [A1 A1*c1+b1;
        A2 A2*c2+b2]
    b = [b1; b2]
    l = -b
    u = fill(Inf, length(l))
    return A, b, P, q
end

# given two polys and two interior points, get point of intersection and signed distance, i.e. a sharing factor s between these two polys
# min s, s.t. A1*(p-c1)+(1+s)(A1*c1+b1)≥0; A2*(p-c2)+(1+s)(A2*c2+b2)≥0 
function get_single_sd_point_shrunk(P1, P2; c1=0, c2=0, eps_abs = 1e-6, eps_rel=1e-6, max_iter=4000)
    A, b, P, q = get_LP_data_point_shrunk(P1, P2; c1=c1, c2=c2)
    l = -b
    u = fill(Inf, length(l))
    #Main.@infiltrate
    ret = solve_qp(UseOSQPSolver(); P=P, q=q, A=A, l=l, u=u, polish=true, verbose=false, eps_abs=eps_abs, eps_rel=eps_rel, max_iter=max_iter)
    #ret.info.status != Symbol("Solved") && @warn ret.info.status
    return ret
end