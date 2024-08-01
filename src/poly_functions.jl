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
                norm(t) < tol && continue
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
    for i in 1:N-1
        vii = @lift($V[i+1])
        vi = @lift($V[i])
        xs = @lift([$vi[1], $vii[1]])
        ys = @lift([$vi[2], $vii[2]])
        lines!(ax, xs, ys; kwargs...)
    end
    vii = @lift($V[1])
    vi = @lift($V[N])
    xs = @lift([$vi[1], $vii[1]])
    ys = @lift([$vi[2], $vii[2]])
    lines!(ax, xs, ys; kwargs...)
end

function get_edge_ind(A, b, vi, vii; tol=1e-4)
    ind_vi = A * vi + b .< tol
    ind_vii = A * vii + b .< tol
    ind_edge = findall(x -> x==2, ind_vi + ind_vii)
    if length(ind_edge) == 0
        return 0
    else
        return ind_edge[1]
    end
end

function delete_first_if_m1(; kwargs...)
    kwargs_list = collect(kwargs)
    m1 = 0
    if length(kwargs_list) > 0 && kwargs_list[1][1] == :m1
        m1 = kwargs_list[1][2]
        deleteat!(kwargs_list, 1)
    end
    new_kwargs = (; kwargs_list...)
    return m1, new_kwargs
end

function plot_with_indices(ax, P::ConvexPolygon2D; kwargs...)
    N = length(P.V)
    V = P.V
    A = P.A
    b = P.b
    m1, kwargs = delete_first_if_m1(; kwargs...)
    for i in 1:N-1
        vii = V[i+1]
        vi = V[i]
        lines!(ax, [vi[1], vii[1]], [vi[2], vii[2]]; kwargs...)

        x_center = (vi[1] + vii[1])/2
        y_center = (vi[2] + vii[2])/2
        ind = get_edge_ind(A, b, vi, vii) + m1
        ind_text = "$ind"
        text!(ax, x_center, y_center; align=(:center, :center), text=ind_text, color=:black, fontsize=20)
    end
    vii = V[1]
    vi = V[N]
    lines!(ax, [vi[1], vii[1]], [vi[2], vii[2]]; kwargs...)

    x_center = (vi[1] + vii[1])/2
    y_center = (vi[2] + vii[2])/2
    ind = get_edge_ind(A, b, vi, vii) + m1
    ind_text = "$ind"
    text!(ax, x_center, y_center; align=(:center, :center), text=ind_text, color=:black, fontsize=20)
end

function plot_with_indices(ax, P::Observable{ConvexPolygon2D}; kwargs...)
    N = length(P[].V)
    V = @lift($P.V)
    A = @lift($P.A)
    b = @lift($P.b)
    m1, kwargs = delete_first_if_m1(; kwargs...)
    
    for i in 1:N-1
        vii = @lift($V[i+1])
        vi = @lift($V[i])
        xs = @lift([$vi[1], $vii[1]])
        ys = @lift([$vi[2], $vii[2]])
        lines!(ax, xs, ys; kwargs...)

        x_center = @lift(($vi[1] + $vii[1])/2)
        y_center = @lift(($vi[2] + $vii[2])/2)
        ind = @lift(get_edge_ind($A, $b, $vi, $vii) + m1)
        ind_text = GLMakie.lift(ind) do in
            "$(in)"
        end
        text!(ax, x_center, y_center; align=(:center, :center), text=ind_text, color=:black, fontsize=20)
    end
    
    vii = @lift($V[1])
    vi = @lift($V[N])
    xs = @lift([$vi[1], $vii[1]])
    ys = @lift([$vi[2], $vii[2]])
    lines!(ax, xs, ys; kwargs...)

    x_center = @lift(($vi[1] + $vii[1])/2)
    y_center = @lift(($vi[2] + $vii[2])/2)
    ind = @lift(get_edge_ind($A, $b, $vi, $vii) + m1)
    ind_text = GLMakie.lift(ind) do in
        "$(in)"
    end
    text!(ax, x_center, y_center; align=(:center, :center), text=ind_text, color=:black, fontsize=20)
end
# PolyPlanning.plot_xt(x0, Pe.A, Pe.b, Po.A, Po.b)
function plot_xt(xt, Ae, be, Ao, bo; fig=Figure(), ax=Axis(fig[1, 1], aspect=DataAspect()))
    xx = xt[1:3]
    Aeb = shift_to(Ae, be, xx)
    self_poly = ConvexPolygon2D(Aeb[1], Aeb[2])
    plot_with_indices(ax, self_poly; color=:blue)
    obs_poly = ConvexPolygon2D(Ao, bo)
    plot_with_indices(ax, obs_poly; m1=length(be), color=:red)
    
    display(GLMakie.Screen(), fig)
    (fig, ax)
end

function plot_xt(xt, ego_polys, obs_polys; fig=Figure(), ax=Axis(fig[1, 1], aspect=DataAspect()))
    xx = xt[1:3]
    for ego_poly in ego_polys
        Aeb = shift_to(ego_poly.A, ego_poly.b, xx)
        self_poly = ConvexPolygon2D(Aeb[1], Aeb[2])
        plot!(ax, self_poly; color=:blue)
    end

    for obs_poly in obs_polys
        plot!(ax, obs_poly; color=:red)
    end
    display(GLMakie.Screen(), fig)
    (fig, ax)
end

function solve_lp(AA, bb, qq)
    ## use JuMP and Clp solver
    model = JuMP.Model(Clp.Optimizer)
    JuMP.set_attribute(model, "LogLevel", 0) # disable printing log
    JuMP.@variable(model, xx[1:3])
    JuMP.@constraint(model, constraint, AA * xx + bb .>= 0)
    JuMP.@objective(model, Min, qq' * xx)
    JuMP.optimize!(model)

    status = JuMP.termination_status(model)
    if status != OPTIMAL
        @warn status
        @infiltrate
    else
        primals = JuMP.value.(xx)
        duals = JuMP.dual.(constraint)
        cons = AA * primals + bb
    end
    primals, duals, cons
end

# PolyPlanning.plot_inflated(x0, Pe.A, Pe.b, Po.A, Po.b)
function plot_inflated(xt, Ae, be, Ao, bo; p_intercept=nothing, sd=nothing, fig=Figure(), ax=Axis(fig[1, 1], aspect=DataAspect()))
    xx = xt[1:3]
    Aex, bex = shift_to(Ae, be, xx)
    self_poly = ConvexPolygon2D(Aex, bex)
    centroidex = sum(self_poly.V)/length(self_poly.V)
    plot!(ax, self_poly; color=:blue)
    obs_poly = ConvexPolygon2D(Ao, bo)
    centroido = sum(obs_poly.V)/length(obs_poly.V)
    plot!(ax, obs_poly; color=:red)
    
    if isnothing(sd)
        AA, bb, qq = gen_LP_data(Aex, bex, centroidex, Ao, bo, centroido)
        primals = solve_lp(AA, bb, qq)[1]
        p_intercept = primals[1:2]
        sd = primals[3]
    end

    # draw inflated
    be_inflated = bex + sd * (Aex * centroidex + bex)
    bo_inflated = bo + sd * (Ao * centroido + bo)
    ego_inflated = ConvexPolygon2D(Aex, be_inflated)
    obs_inflated = ConvexPolygon2D(Ao, bo_inflated)

    plot_with_indices(ax, ego_inflated; color=:lightblue, linewidth=2)
    plot_with_indices(ax, obs_inflated; m1=length(be), color=:pink, linewidth=2)

    # draw the intercept and centroids
    scatter!(ax, centroidex[1], centroidex[2]; color=:blue)
    scatter!(ax, centroido[1], centroido[2]; color=:red)
    if !isnothing(p_intercept)
        scatter!(ax, p_intercept[1], p_intercept[2]; color=:green)
    end

    display(GLMakie.Screen(), fig)
    (fig, ax)
end

# check if a halfspace ai'x+bi≥0 already exists in Ax+b≥0
function if_halfspace_exist(A, b, ai, bi; tol=1e-5)
    if length(A) == 0
        return false
    end
    shapeA = (length(A[1]), length(A))
    matrixA = reshape(vcat(A...), shapeA)'
    err = hcat(matrixA, b) .- [ai; bi]'
    row_nms = [norm(row) for row in eachrow(err)]
    if sum(row_nms .< tol) > 0
        return true
    end
    false
end

# check if a vertex already exists in V
function if_vertex_exist(V, v; tol=1e-5)
    if length(V) == 0
        return false
    end
    shapeV = (length(V[1]), length(V))
    matrixV = reshape(vcat(V...), shapeV)'
    err = matrixV .- v'
    row_nms = [norm(row) for row in eachrow(err)]
    if sum(row_nms .< tol) > 0
        return true
    end
    false
end

struct ConvexPolygon3D
    A::SparseMatrixCSC{Float64, Int64}
    b::Vector{Float64}
    V::Vector{Vector{Float64}}
    c::Vector{Float64}
    function ConvexPolygon3D(A,b,V,c)
        new(A,b,V,c)
    end
    function ConvexPolygon3D(A,b; tol=1e-5)
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
                for i3 = i2+1:m
                    a3 = A[i3,:]
                    rhs = [b[i1], b[i2], b[i3]]
                    try
                        v = [a1'; a2'; a3'] \ -rhs
                        if all(- tol .≤ A*v+b) && !if_vertex_exist(V, v)
                            push!(V, v)
                        end
                    catch
                    end
                end
            end
        end
        n_verts = length(V)
        c = sum(V) ./ n_verts
    
        θs = map(V) do vi
            d = vi-c
            atan(d[3], sqrt(d[1]^2 + d[2]^2))
        end
        I = sortperm(θs) |> reverse
        V = V[I]

        new(A,b,V,c)
    end
    function ConvexPolygon3D(V; tol=1e-5)
        # Convert from Vrep to Hrep
        A = []
        b = Float64[]
        N = length(V)
        supporting_verts = Int[]
        for i = 1:N
            v1 = V[i]
            for j = i+1:N
                v2 = V[j]
                t1 = v2 - v1
                norm(t1) < tol && continue
                for k = j+1:N
                    v3 = V[k]
                    t2 = v3 - v1
                    norm(t2) < tol && continue
                    t1 ./= norm(t1)
                    t2 ./= norm(t2)
                    ai = cross(t1, t2)
                    norm(ai) ≤ tol && continue
                    ai ./= norm(ai)
                    bi = -(ai'*v1)
                    if all( -tol .≤ ai'*v+bi for v in V)
                        push!(supporting_verts, i)
                        push!(supporting_verts, j)
                        push!(supporting_verts, k)
                        if !if_halfspace_exist(A, b, ai, bi)
                            push!(A, ai)
                            push!(b, bi)
                        end
                    elseif all( tol .≥ ai'*v+bi for v in V)
                        push!(supporting_verts, i)
                        push!(supporting_verts, j)
                        push!(supporting_verts, k)
                        if !if_halfspace_exist(A, b, -ai, -bi)
                            push!(A, -ai)
                            push!(b, -bi)
                        end
                    end
                end
            end
        end
        
        A = hcat(A...)' 
        θs = map(1:size(A,1)) do i
            ai = A[i,:]
            θ = atan(ai[3], sqrt(ai[1]^2 + ai[2]^2))
        end
        I = sortperm(θs) |> reverse 
        A = A[I,:]
        b = b[I]
        
        supporting_verts = Set(supporting_verts) |> collect
        V = V[supporting_verts]
        c = sum(V) / length(V)
        θs = map(V) do vi
            d = vi-c
            atan(d[3], sqrt(d[1]^2 + d[2]^2))
        end
        I = sortperm(θs) |> reverse
        V = V[I]
        new(sparse(A),b,V,c)
    end
end

function getAb_from_polyhedra(poly)
    polyh = hrep(poly)
    A = []
    b = Float64[]
    for i in eachindex(polyh.halfspaces)
        hss = polyh.halfspaces[i]
        ai = -hss.a
        nm = norm(ai)
        ai ./= nm
        bi = hss.β / nm
        push!(A, ai)
        push!(b, bi)
    end
    A = hcat(A...)'
    return sparse(A), b
end

function getV_from_polyhedra(poly)
    polyv = vrep(poly)
    return polyv.points.points
end

function get_polyhedra_from_V(V)
    polyv = vrep(V)
    polyhedron(polyv)
end


function get_polyhedra_from_Ab(A, b)
    hss = map(1:length(b)) do i
        HalfSpace(-A[i, :], b[i])
    end
    hr = hrep(hss)
    polyhedron(hr)
end

function plot_3D!(ax3, P::ConvexPolygon3D; kwargs...)
    poly = get_polyhedra_from_V(P.V)
    mesh_P = Polyhedra.Mesh(poly)
    GLMakie.wireframe!(ax3.scene, mesh_P; kwargs...)
end

function plot_3D!(ax3, P::Observable{ConvexPolygon3D}; kwargs...)
    poly = @lift(get_polyhedra_from_V($P.V))
    mesh_P = @lift(Polyhedra.Mesh($poly))
    GLMakie.wireframe!(ax3.scene, mesh_P; kwargs...)
end

# function plot_3D(ax3, P::Observable{ConvexPolygon3D}; kwargs...)
#     mesh_P = @lift(Polyhedra.Mesh($P))
#     GLMakie.wireframe!(ax3.scene, mesh_P; kwargs...)
# end

# if three vertices are on the same face, return face id; otherwise, return 0
function get_face_ind(A, b, vi, vii, viii; tol=1e-4)
    ind_vi = A * vi + b .< tol
    ind_vii = A * vii + b .< tol
    ind_viii = A * viii + b .< tol
    ind_face = findall(x -> x==3, ind_vi + ind_vii + ind_viii)
    if length(ind_face) == 0
        return 0
    else
        return ind_face[1]
    end
end

function plot_with_indices(ax3, P::ConvexPolygon3D; kwargs...)
    N = length(P.V)
    V = P.V
    A = P.A
    b = P.b
    m1, kwargs = delete_first_if_m1(; kwargs...)
    plot_3D!(ax3, P; kwargs...)

    for i = 1:N
        vi = V[i]
        for j = i+1:N
            vii = V[j]
            for k = j+1:N
                viii = V[k]
                ind = get_face_ind(A, b, vi, vii, viii) + m1
                if !(ind==m1)
                    ind_text = "$ind"
                    x_center = (vi[1] + vii[1]+ viii[1])/3
                    y_center = (vi[2] + vii[2]+ viii[2])/3
                    z_center = (vi[3] + vii[3]+ viii[3])/3
                    text!(ax3, x_center, y_center, z_center; align=(:center, :center), text=ind_text, color=:black, fontsize=20)
                end
            end
        end
    end
end

function mrp_from_axis_angle(e, θ)
    return tan(θ/4)*e
end

function axis_angle_from_mrp(mrp)
    nm = norm(mrp)
    θ = 4 * atan(nm)
    return mrp/nm, θ
end

function mrp_from_q(q)
    return q[2:4] / (1+q[1])
end

function q_from_mrp(mrp)
    nm2 = mrp'mrp
    return (1/(1+nm2))*[1-nm2; 2*mrp]
end

function R_from_q(q; tol=1e-4)
    if abs(norm(q)-1)> tol
        @warn "not a unit quaternion"
        @infiltrate
        q = normalize(q)
    end
    q1, q2, q3, q4 = q
    
    [(2*q2^2+2*q1^2-1)   2*(q2*q3 - q4*q1)   2*(q2*q4 + q3*q1)
    2*(q2*q3 + q4*q1)  (2*q3^2+2*q1^2-1)   2*(q3*q4 - q2*q1)
    2*(q2*q4 - q3*q1)   2*(q3*q4 + q2*q1)  (2*q4^2+2*q1^2-1)]
end

function skew(v)
    [0 -v[3] v[2]
    v[3] 0 -v[1]
    -v[2] v[1] 0]
end

function R_from_mrp(mrp)
    nm2 = mrp'*mrp
    I + ( 4*(1-nm2)*I + 8*skew(mrp) ) * skew(mrp) / (1+nm2)^2
end

# x∈R12, x[1:3] is position, x[4:6] is orientation, x[7:9] is velocity, x[10:12] is angular velocity
function shift_to_3D(A, b, x)
    p = x[1:3]
    mrp = x[4:6]
    R = R_from_mrp(mrp)
    At = A * R'
    bt = b - At * p
    At, bt
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
