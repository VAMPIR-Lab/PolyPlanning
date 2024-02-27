struct ConvexPolygon2D
    A::SparseMatrixCSC{Float64, Int64}
    b::Vector{Float64}
    V::Vector{Vector{Float64}}
    c::Vector{Float64}
    function ConvexPolygon2D(A,b; tol=1e-5)
        # Convert from Hrep to Vrep
        m = length(b)
        V = Vector{Float64}[]
        for i in 1:m
            A[i,:] ./= norm(A[i,:])
            b[i] ./= norm(A[i,:])
        end
        for i1 = 1:m
            a1 = A[i1,:]
            for i2 = i1+1:m
                a2 = A[i2,:]
                rhs = [b[i1], b[i2]]
                try
                    v = [a1'; a2'] \ rhs
                    if all(- tol .≤ A*v+b)
                        push!(V, v)
                    end
                catch
                end
            end
        end
        n_verts = length(V)
        c = sum(V) ./ n_verts
    
        θs = map(V) do vi
            d = vi-c
            atan(d[2],d[1])
        end
        I = sortperm(θs) |> reverse
        V = V[I]

        new(A,b,V,c)
    end
    function ConvexPolygon2D(V; tol=1e-5)
        # Convert from Vrep to Hrep
        c = sum(V) / length(V)
        θs = map(V) do vi
            d = vi-c
            atan(d[2],d[1])
        end
        I = sortperm(θs) |> reverse
        V = V[I]

        A = []
        b = Float64[]
        N = length(V)
        for i = 1:N-1
            v1 = V[i]
            v2 = V[i+1]
            t = v2-v1
            t ./= norm(t)
            ai = [t[2], -t[1]]
            bi = -(ai'*v1)
            push!(A,ai)
            push!(b,bi)
        end
        v1 = V[N]
        v2 = V[1]
        t = v2-v1
        t ./= norm(t)
        ai = [t[2], -t[1]]
        bi = -(ai'*v1)
        push!(A,ai)
        push!(b,bi)
        A = hcat(A...)' 
        c = sum(V) ./ N
        new(sparse(A),b,V,c)
    end
end

function plot!(ax, P::ConvexPolygon2D; kwargs...)
    N = length(P.V) 
    V = P.V
    c = P.c
    for i in 1:N-1
        vii = V[i+1]
        vi = V[i]
        lines!(ax, [vi[1], vii[1]], [vi[2], vii[2]]; kwargs...)
    end
    vii = V[1]
    vi = V[N]
    lines!(ax, [vi[1], vii[1]], [vi[2], vii[2]]; kwargs...)
end

function signed_distance(P1::ConvexPolygon2D, 
                         P2::ConvexPolygon2D)
    A1 = -P1.A
    b1 = -P1.b
    A2 = P2.A
    b2 = P2.b

    m1 = length(b1)
    m2 = length(b2)
    M = [spzeros(m1,m1) A1 spzeros(m1,m2) sparse(-ones(m1,1));
         -A1'  spzeros(2, 2) -A2' spzeros(2,1);
         spzeros(m2, m1) A2 spzeros(m2,m2) spzeros(m2,1);
         sparse(ones(1,m1)) spzeros(1,2) spzeros(1,m2) spzeros(1,1)]
    q = [b1; zeros(2); b2; -1]
    l = [zeros(m1); fill(-Inf, 2); zeros(m2); fill(-Inf, 1)]
    u = fill(Inf, m1+m2+3)
    ret = solve_lmcp(UsePATHSolver(), M, q, l, u) 
    y = ret.z[1:m1]
    x = ret.z[m1+1:m1+2]
    sd = y'*(A1*x+b1)
    (; sd, y, x, ret.solve_status, ret.info)
end
