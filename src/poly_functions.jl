struct ConvexPolygon2D
    A::SparseMatrixCSC{Float64, Int64}
    l::Vector{Float64}
    u::Vector{Float64}
    V::Vector{Vector{Float64}}
    c::Vector{Float64}
    function ConvexPolygon2D(A,l,u; tol=1e-5)
        # Convert from Hrep to Vrep
        m = length(l)
        V = Vector{Float64}[]
        for i1 = 1:m
            a1 = A[i1,:]
            for i2 = i1+1:m
                a2 = A[i2,:]
                rhs = [l[i1] l[i1] u[i1] u[i1];
                       l[i2] u[i2] l[i2] u[i2]]
                try
                    v = [a1'; a2'] \ rhs
                    for j in 1:4
                        vj = v[:,j]
                        if all(l .- tol .≤ A*vj .≤ u .+ tol)
                            push!(V, vj)
                        end
                    end
                catch
                end
            end
        end
        n_verts = length(V)
        c = sum(V) ./ n_verts
        new(A,l,u,V,c)
    end
    function ConvexPolygon2D(V; tol=1e-5)
        # Convert from Vrep to Hrep
        A = []
        l = Float64[]
        u = Float64[]
        N = length(V)
        for i1 = 1:N
            v1 = V[i1]
            for i2 = i1+1:N
                v2 = V[i2]
                t = v2-v1
                ai = [-t[2], t[1]]
                li = ai'*v1
                if all(ai'*v ≥ li - tol for v in V)
                    push!(A,ai)
                    push!(l,li)
                else
                    ai = -ai
                    li = ai'*v1
                    if all(ai'*v ≥ li - tol for v in V)
                        push!(A,ai)
                        push!(l,li)
                    end
                end
            end
        end
        A = hcat(A...)' 
        u = fill(Inf, length(l))
        c = sum(V) ./ N
        new(sparse(A),l,u,V,c)
    end
end

function plot!(ax, P::ConvexPolygon2D; kwargs...)
    N = length(P.V) 
    V = P.V
    c = P.c
    VV = Vector{Float64}[]
    θs = map(V) do vi
        d = vi-c
        atan(d[2],d[1])
    end
    I = sortperm(θs) |> reverse
    V = V[I]
    for i in 1:N-1
        vii = V[i+1]
        vi = V[i]
        lines!(ax, [vi[1], vii[1]], [vi[2], vii[2]]; kwargs...)
    end
    vii = V[1]
    vi = V[N]
    lines!(ax, [vi[1], vii[1]], [vi[2], vii[2]]; kwargs...)
    scatter!(ax, [P.c[1],], [P.c[2],]; kwargs...)
end

function signed_distance(P1::ConvexPolygon2D, 
                         P2::ConvexPolygon2D;
                         ax = nothing,
                         phase_1_tol=1e-6)
    I2 = sparse(1.0I, 2, 2)
    P = [I2 -I2;
            -I2 I2]
    m1 = length(P1.l)
    m2 = length(P2.l)
    Z1 = spzeros(m1,2)
    Z2 = spzeros(m2,2)
    A = [P1.A Z1;
            Z2 P2.A]
    l = [P1.l; P2.l]
    u = [P1.u; P2.u]
    ret = solve_qp(OSQPSolver(); P, A, l, u, verbose=false, polish=true)
    
    if ret.info.obj_val > phase_1_tol
        if !isnothing(ax)
            x1 = ret.x[1:2]
            x2 = ret.x[3:4]
            scatter!(ax, [x1[1],], [x1[2],], color=:black)
            scatter!(ax, [x2[1],], [x2[2],], color=:black)
        end
        return (; sd = ret.info.obj_val, ret)
    else
        
        V = P2.V
        c = P2.c
        VV = Vector{Float64}[]
        θs = map(V) do vi
            d = vi-c
            atan(d[2],d[1])
        end
        I = reverse(sortperm(θs))
        V = V[I]
        N = length(I)
        for i in 1:N-1
            t = V[i+1]-V[i]
            t ./= norm(t)
            n = [-t[2], t[1]]
        end
        return (; sd=0.0)
        ## phase 2
        #d = P1.c-P2.c
        #d ./= norm(d)
        #q = [d; -d] 
        #A = [P1.A Z1;
        #     P2.A Z2;
        #     Z1 P1.A;
        #     Z2 P2.A]
        #l = [l; l]
        #u = [u; u]
        #ret = solve_qp(OSQPSolver(); q, A, l, u, verbose=false, polish=true, eps_abs=1e-8, eps_rel=1e-8)
        #if !isnothing(ax)
        #    x1 = ret.x[1:2]
        #    x2 = ret.x[3:4]
        #    scatter!(ax, [x1[1],], [x1[2],], color=:black)
        #    scatter!(ax, [x2[1],], [x2[2],], color=:black)
        #end
        #return (; sd = ret.info.obj_val, ret)
    end
end

function signed_distance(P::ConvexPolygon2D, p::Vector{Float64})
    V = P.V
    c = P.c
    VV = Vector{Float64}[]
    θs = map(V) do vi
        d = vi-c
        atan(d[2],d[1])
    end
    I = reverse(sortperm(θs))
    V = V[I]
    N = length(I)
    sds = Float64[]
    for i in 1:N-1
        t = V[i+1]-V[i]
        t ./= norm(t)
        ni = [-t[2], t[1]]
        push!(sds, ni'*(p-V[i]))
    end
    t = V[1]-V[N]
    t ./= norm(t)
    ni = [-t[2], t[1]]
    push!(sds, ni'*(p-V[N]))
    id = argmax(sds)
    return (; sd=sds[id], id)
end

function test_psd()
    P = ConvexPolygon2D([[0.0,0], [1,0], [-0.5,1], [1.5,1], [0.5,2.5]])
    G = [Vector{Float64}[] for _ in 1:5] 
    H = [Vector{Float64}[] for _ in 1:5] 
    res = 0.01
    for t1 in -1:res:2
        for t2 in -0.5:res:3
            t = [t1,t2]
            ret = signed_distance(P, t)
            if ret.sd > 0
                push!(G[ret.id], t)
            else
                push!(H[ret.id], t)
            end
        end
    end
    f = Figure()
    ax = Axis(f[1,1], aspect=DataAspect())
    colors = [:blue, :red, :green, :yellow, :pink]
    for i in 1:5
        xpts = [g[1] for g in G[i]]
        ypts = [g[2] for g in G[i]]
        scatter!(ax, xpts, ypts, color=colors[i], markersize=2)
        xpts = [g[1] for g in H[i]]
        ypts = [g[2] for g in H[i]]
        scatter!(ax, xpts, ypts, color=colors[i], marker='*', markersize=2)
    end
    plot!(ax, P; color=:black)

    display(f)
end


function test_sd(;trange=-5:0.1:5, plot_iters=[])
    P1 = ConvexPolygon2D([[1.0,0],[0,2],[-1,0],[0,-1]])
    sds = Float64[]
    for (e,t) in enumerate(trange)
        P2 = ConvexPolygon2D([[-0.2,3+t],[-0.1,-0.5+t],[1.3,2+t]])
        if e in plot_iters
            f = Figure(); ax = Axis(f[1,1], aspect=DataAspect())
            plot!(ax, P1; color=:blue)
            plot!(ax, P2; color=:red)
            lines!(ax, [P1.c[1], P2.c[1]], [P1.c[2], P2.c[2]], color=:black)
            (; sd) = signed_distance(P1, P2; ax)
            display(f)
        else
            (; sd) = signed_distance(P1, P2)
        end
        push!(sds, sd)
    end
    f = Figure(); ax = Axis(f[1,1])
    scatter!(ax, trange, sds)
    display(f)
end
