using Revise
using PolyPlanning
using LinearAlgebra
using OSQP
using SparseArrays
using GLMakie
using PATHSolver
using Symbolics
using ProgressMeter
using Random
using Combinatorics
using ArgCheck
using Statistics
using JLD2
using JuMP
using Clp
using Polyhedra
using Infiltrator

function generate_s_t_plot(P, O, c, x1, x2; num_interp = 100, interval_display = 10, eps_abs=1e-6, eps_rel=1e-6, max_iter=4000)

    # plot the obstacle
    fig = PolyPlanning.Figure()
    ax = PolyPlanning.Axis(fig[1,1], aspect=PolyPlanning.DataAspect())
    PolyPlanning.plot!(ax, O)

    # plot the starting poly and the ending poly
    P1 = PolyPlanning.shift_to_around_point(P, x1; c=c)
    P2 = PolyPlanning.shift_to_around_point(P, x2; c=c)
    PolyPlanning.plot!(ax, P1)
    PolyPlanning.plot!(ax, P2)

    # plot interpolated polys
    polys = PolyPlanning.interpolation_polys(P, x1, x2; c=c, num = num_interp)
    
    sds = [] # signed distance array
    ps = [] # position of intersection point array
    for i in 1:num_interp
        if i % interval_display == 1 || i == num_interp
            PolyPlanning.plot!(ax, polys[i])
        end
        ret = PolyPlanning.get_single_sd_point_shrunk(polys[i], O; eps_abs=eps_abs, eps_rel=eps_rel, max_iter=max_iter)
        ret.info.status != Symbol("Solved") && @warn "not solved accurately" sattus_type=ret.info.status index_interp=i
        append!(sds, ret.x[end])
        append!(ps, [ret.x[1:end-1]])
    end

    t = [i/(num_interp-1) for i in 0:num_interp-1]
    # plot s-t curve
    fig2 = Figure()
    ax2 = Axis(fig2[1, 1], title = "s-t Plot", xlabel = "t", ylabel = "s")
    scatter!(ax2, t, sds, markersize = 5, color = :red)

    # plot trajectory of intersection point p
    fig3 = Figure()
    ax3 = Axis(fig3[1, 1], title = "trajectory of intersection point")#, xlabel = "t", ylabel = "s")
    px = [point[1] for point in ps]
    py = [point[2] for point in ps]
    scatter!(ax3, px, py, markersize = 5, color = :red)
    PolyPlanning.plot!(ax3, O)
    PolyPlanning.plot!(ax3, P1)
    PolyPlanning.plot!(ax3, P2)

    # create new windows to display plots
    display(GLMakie.Screen(), fig)
    display(GLMakie.Screen(), fig2)
    display(GLMakie.Screen(), fig3)
end


function plot_sd_point_shrunk(P1, P2; c1=0, c2=0, eps_abs=1e-6, eps_rel=1e-6, max_iter=4000)
    ret = PolyPlanning.get_single_sd_point_shrunk(P1, P2; c1=c1, c2=c2, eps_abs=eps_abs, eps_rel=eps_rel, max_iter=max_iter)
    p = ret.x[1:end-1]
    s = ret.x[end]
    @info "intersection point" p
    @info "signed distance" s
    ret.info.status != Symbol("Solved") && @warn "not solved accurately" indsattus_type=ret.info.status
    P1_shrunk = PolyPlanning.ConvexPolygon2DPointShrunk(P1; c=c1, s)
    P2_shrunk = PolyPlanning.ConvexPolygon2DPointShrunk(P2; c=c2, s)

    fig = PolyPlanning.Figure()
    ax = PolyPlanning.Axis(fig[1,1], aspect=PolyPlanning.DataAspect())
    fig2 = PolyPlanning.Figure()
    ax2 = PolyPlanning.Axis(fig2[1,1], aspect=PolyPlanning.DataAspect())
    PolyPlanning.plot!(ax, P1)
    PolyPlanning.plot!(ax, P2)
    PolyPlanning.plot!(ax2, P1_shrunk)
    PolyPlanning.plot!(ax2, P2_shrunk)
    display(GLMakie.Screen(), fig)
    display(GLMakie.Screen(), fig2)
    return ret
end

function ConvexPolygon2DPointShrunk(P; c=0, s=0)
    A = P.A
    b = P.b
    V = P.V
    if c==0
        c = sum(V) / length(V)
    end
    b = b + s * (A * c + b)
    PolyPlanning.ConvexPolygon2D(A, b)
end

# # generate polygons given a series of factor s; when s = -1, the polygon shrinks to a point, which may cause some error when ploting
# function GenerateShrinkingPolys(P; c=0, S=[-0.99, -0.9, -0.8, -0.6, -0.3, 0])
#     ShrunkPolys = [ConvexPolygon2DPointShrunk(P; c, s) for s in S]
#     return ShrunkPolys
# end

function gen_LP_data(xt, Ae::AbstractArray{T}, be, centroide, Ao, bo, centroido; is_newsd=false) where {T}
    A = [Ae Ae*centroide+be
        Ao Ao*centroido+bo]
    b = [be; bo]
    q = [0, 0, 1.0]
    # print("new")
    (A, b, q)
end

function if_ass_feasible(ass, m1, m2)
    if_ego = false
    if_obs = false
    for i in ass
        if i ∈ [i for i in 1:m1]
            if_ego=true
        end
        if i ∈ [i for i in m1+1:m1+m2]
            if_obs=true
        end
    end
    return if_ego && if_obs
end

function g_col_single(xt, Ae, be, centroide, Ao, bo, centroido; is_newsd=false)
    sds = Dict()
    Aex, bex = PolyPlanning.shift_to(Ae, be, xt)
    R = [cos(xt[3]) sin(xt[3])
        -sin(xt[3]) cos(xt[3])]
    centroidex = xt[1:2] + R * centroide
    AA, bb, qq = gen_LP_data(xt, Aex, bex, centroidex, Ao, bo, centroido; is_newsd=is_newsd)
    m1 = length(bex)
    m2 = length(bo)
    M = [zeros(Num, 3, 3) -AA'
        AA zeros(Num, m1 + m2, m1 + m2)]
    r = [qq; bb]
    all_active_inds = collect(1:m1+m2)
    Itr = powerset(all_active_inds) |> collect
    Itr_reduced=[]
    for ass in Itr
        if if_ass_feasible(ass, m1, m2)
            push!(Itr_reduced, ass)
        end
    end
    Itr = copy(Itr_reduced)

    for active_inds in Itr
        length(active_inds) > 3 && continue
        assignment = [i ∈ active_inds for i in 1:m1+m2]
        try
            AA_active = collect(AA[active_inds, :])
            bb_active = collect(bb[active_inds])
            # TODO what if not unique primal? Need to resolve
            if length(active_inds) == 3
                zz = -AA_active \ bb_active
                #zz = -(AA_active'*AA_active)\AA_active'*bb_active
            else
                #if linear system is underdetermined, use minimum norm solution (calculated by right inverse)
                #Note: every solution has the same sd, i.e., zz[3], but different zz[1:2]
                zz = -AA_active' * ((AA_active * AA_active') \ bb_active)
            end
            sd = zz#[3]
            sds[active_inds] = sd
        catch err
            if err isa LinearAlgebra.SingularException
                continue
            else
                @warn(err)
            end
        end
    end
    sds, AA, bb, qq
end

function get_single_sd_ids(xt, Ae, be, centroide, Ao, bo, centroido, max_derivs, t; is_newsd=false)
    Aex, bex = PolyPlanning.shift_to(Ae, be, xt)
    R = [cos(xt[3]) sin(xt[3])
        -sin(xt[3]) cos(xt[3])]
    centroidex = xt[1:2] + R * centroide
    AA, bb, qq = gen_LP_data(xt, Aex, bex, centroidex, Ao, bo, centroido; is_newsd=is_newsd)
    m1 = length(bex)
    m2 = length(bo)

    ## use JuMP and Clp solver
    model = JuMP.Model(Clp.Optimizer)
    JuMP.set_attribute(model, "LogLevel", 0) # disable printing log
    JuMP.@variable(model, xx[1:3])
    JuMP.@constraint(model, constraint, AA * xx + bb .>= 0)
    JuMP.@objective(model, Min, qq' * xx)
    JuMP.optimize!(model)

    status = JuMP.termination_status(model)
    if status != OPTIMAL
        duals = zeros(m1 + m2)
        cons = duals
        xx = JuMP.value.(xx)
        @warn status
        #Main.@infiltrate
        #plot_inflated(xt, Ae, be, centroidex, Ao, bo, centroido, xx[1:2], xx[3])
    else
        primals = JuMP.value.(xx)
        duals = JuMP.dual.(constraint)
        cons = AA * primals + bb
    end

    tol = 1e-3
    I1 = duals .≥ tol .&& cons .< tol
    I2 = duals .< tol .&& cons .< tol
    I3 = duals .< tol .&& cons .≥ tol

    return I1, I2, I3, primals, duals, cons
end


Ve = [[-1, -1], [-1.1, 1.2], [0.6, 0.8], [1.7, -0.5]]
Pe = PolyPlanning.ConvexPolygon2D(Ve)
centroide = sum(Pe.V) / length(Pe.V)
Ae = Pe.A
be = Pe.b

Vo = [[-1, -1], [-1, 1.0], [1, 1], [1, -1]]
Po = PolyPlanning.ConvexPolygon2D(Ve)
centroido = sum(Po.V) / length(Po.V)
Ao = Po.A
bo = Po.b

xt=[2, 2, π*0, 0, 0, 0]
Aex, bex = PolyPlanning.shift_to(Ae, be, xt)
Pex = PolyPlanning.ConvexPolygon2D(Aex, bex)
centroidex = sum(Pex.V) / length(Pex.V)


fig = PolyPlanning.Figure()
ax = PolyPlanning.Axis(fig[1,1], aspect=PolyPlanning.DataAspect())
PolyPlanning.plot!(ax, Pex; color=:blue)
PolyPlanning.plot!(ax, Po; color=:red)



sds, AA, bb, qq = g_col_single(xt, Ae, be, centroide, Ao, bo, centroido; is_newsd=false)
I1, I2, I3, primals, duals, cons = get_single_sd_ids(xt, Ae, be, centroide, Ao, bo, centroido, 4, 20; is_newsd=false)

println("sds≥0 and corresponding assignments")
for (i, val) in enumerate(sds)
    if val[2][3]≥-1
        println(val[1], " ",val[2][3])
    end
end


function plot_shrk(Pex, Po, s)
    Pex_shrk = ConvexPolygon2DPointShrunk(Pex; c=0, s=s)
    Po_shrk = ConvexPolygon2DPointShrunk(Po; c=0, s=s)
    PolyPlanning.plot!(ax, Pex_shrk; color=:lightblue)
    PolyPlanning.plot!(ax, Po_shrk; color=:pink)
    display(GLMakie.Screen(), fig)
end

sds_k = collect(keys(sds))
sds_val = collect(values(sds))
feasible_indices = findall(x -> x[3] >= -1, sds_val)
# i = feasible_indices[11]
# println("ass = ", sds_k[i], "\nsigned distance = ", sds_val[i])
# plot_shrk(Pex, Po, sds_val[i][3])

println("\nassignment and corresponding sd which satisifies all 8 constraints:")
tol = 1e-4
for val in sds_val
    err = AA * val + bb 
    if sum(err .>= -tol) == 8
        println("\nprimal vector = ", val)
        println("value of constraints = ", err)
    end
end


hss = map(1:8) do i
    HalfSpace(-AA[i, :], bb[i])
end
hr = hrep(hss) # H-representation for polyhedron 
poly = Polyhedra.polyhedron(hr)
vr = vrep(poly) # V-representation


mesh_poly = Polyhedra.Mesh(poly)
Makie.mesh(mesh_poly, color=:blue)
Makie.wireframe(mesh_poly)