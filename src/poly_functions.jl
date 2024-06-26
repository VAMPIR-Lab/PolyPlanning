struct ConvexPolygon2D
    A::SparseMatrixCSC{Float64, Int64}
    b::Vector{Float64}
    V::Vector{Vector{Float64}}
    function ConvexPolygon2D(A,b,V)
        new(A,b,V)
    end
    function ConvexPolygon2D(A,b; tol=1e-5)
        # Convert from Hrep to Vrep
        m = length(b)
        V = Vector{Float64}[]
        for i in 1:m
            nm = norm(A[i,:])
            A[i,:] ./= nm
            b[i] /= nm
        end
        for i1 = 1:m
            a1 = A[i1,:]
            for i2 = i1+1:m
                a2 = A[i2,:]
                rhs = [b[i1], b[i2]]
                try
                    v = [a1'; a2'] \ -rhs
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

        new(A,b,V)
    end
    function ConvexPolygon2D(V; tol=1e-5)
        # Convert from Vrep to Hrep

        A = []
        b = Float64[]
        N = length(V)
        supporting_verts = Int[]
        for i = 1:N
            for j = 1:N
                j == i && continue
                v1 = V[i]
                v2 = V[j]
                t = v2-v1
                t ./= norm(t)
                ai = [t[2], -t[1]]
                bi = -(ai'*v1)
                if all( -tol .≤ ai'*v+bi for v in V)
                    push!(A,ai)
                    push!(b,bi)
                    push!(supporting_verts, i)
                    push!(supporting_verts, j)
                end
            end
        end
        
        A = hcat(A...)' 
        θs = map(1:size(A,1)) do i
            ai = A[i,:]
            θ = atan(ai[2],ai[1])
        end
        I = sortperm(θs) |> reverse 
        A = A[I,:]
        b = b[I]
        
        supporting_verts = Set(supporting_verts) |> collect
        V = V[supporting_verts]
        c = sum(V) / length(V)
        θs = map(V) do vi
            d = vi-c
            atan(d[2],d[1])
        end
        I = sortperm(θs) |> reverse
        V = V[I]
        new(sparse(A),b,V)
    end
end

function plot!(ax, P::ConvexPolygon2D; kwargs...)
    N = length(P.V) 
    V = P.V
    for i in 1:N-1
        vii = V[i+1]
        vi = V[i]
        lines!(ax, [vi[1], vii[1]], [vi[2], vii[2]]; kwargs...)
    end
    vii = V[1]
    vi = V[N]
    lines!(ax, [vi[1], vii[1]], [vi[2], vii[2]]; kwargs...)
end
function plot!(ax, P::Observable{ConvexPolygon2D}; kwargs...)
    N = length(P[].V)
    V = @lift($P.V)
    A = @lift($P.A)
    b = @lift($P.b)
    function get_edge_ind(A, b, vi, vii)
        ind_vi = A * vi + b .< 1e-4
        ind_vii = A * vii + b .< 1e-4
        ind_edge = findall(x -> x==2, ind_vi + ind_vii)
        if length(ind_edge) == 0
            return 0
        else
            return ind_edge[1]
        end
    end
    
    for i in 1:N-1
        vii = @lift($V[i+1])
        vi = @lift($V[i])
        xs = @lift([$vi[1], $vii[1]])
        ys = @lift([$vi[2], $vii[2]])
        x_center = @lift(($vi[1] + $vii[1])/2)
        y_center = @lift(($vi[2] + $vii[2])/2)
        ind = @lift(get_edge_ind($A, $b, $vi, $vii))
        ind_text = GLMakie.lift(ind) do in
            "$(in)"
        end
        lines!(ax, xs, ys; kwargs...)
        text!(ax, x_center, y_center; align=(:center, :center), text=ind_text, color=:black, fontsize=10)
    end
    
    vii = @lift($V[1])
    vi = @lift($V[N])
    xs = @lift([$vi[1], $vii[1]])
    ys = @lift([$vi[2], $vii[2]])
    x_center = @lift(($vi[1] + $vii[1])/2)
    y_center = @lift(($vi[2] + $vii[2])/2)
    ind = @lift(get_edge_ind($A, $b, $vi, $vii))
    ind_text = GLMakie.lift(ind) do in
        "$(in)"
    end
    lines!(ax, xs, ys; kwargs...)
    text!(ax, x_center, y_center; align=(:center, :center), text=ind_text, color=:black, fontsize=10)
end

function signed_distance(P1::ConvexPolygon2D, 
                         P2::ConvexPolygon2D)
    A1 = P1.A
    b1 = P1.b
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
    z = ret.z
    num_weak = sum(z .< (l .+ 1e-4) .&& (M*z+q .< 1e-4))
    y = ret.z[1:m1]
    x = ret.z[m1+1:m1+2]
    sd = y'*(A1*x+b1)
    (; sd, y, x, ret.solve_status, ret.info, num_weak)
end

function sym_signed_distance(P1::ConvexPolygon2D,
                             P2::ConvexPolygon2D)
    signed_distance(P1,P2) + signed_distance(P2,P1)
end
