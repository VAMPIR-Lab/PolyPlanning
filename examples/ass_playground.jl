using PolyPlanning
using GLMakie

function get_single_sd_ids(xt, Ae, be, centroide, Ao, bo, centroido)
    Aex, bex = PolyPlanning.shift_to(Ae, be, xt)
    R = [cos(xt[3]) sin(xt[3])
        -sin(xt[3]) cos(xt[3])]
    centroidex = xt[1:2] + R * centroide
    AA, bb, qq = PolyPlanning.gen_LP_data(xt, Aex, bex, centroidex, Ao, bo, centroido)
    m1 = length(bex)
    m2 = length(bo)

    # use JuMP and Clp solver
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
    else
        primals = JuMP.value.(xx)
        duals = JuMP.dual.(constraint)
        cons = AA * primals + bb
    end

    tol = 1e-3
    I1 = duals .≥ tol .&& cons .< tol
    I2 = duals .< tol .&& cons .< tol
    I3 = duals .< tol .&& cons .≥ tol

    return primals, duals, cons, I1, I2, I3
end


#function g_col_single(xt, Ae, be, centroide, Ao, bo, centroido; is_newsd=false)
#    sds = Dict()
#    Aex, bex = PolyPlanning.shift_to(Ae, be, xt)
#    R = [cos(xt[3]) sin(xt[3])
#        -sin(xt[3]) cos(xt[3])]
#    centroidex = xt[1:2] + R * centroide
#    AA, bb, qq = PolyPlanning.gen_LP_data(xt, Aex, bex, centroidex, Ao, bo, centroido; is_newsd=is_newsd)
#    m1 = length(bex)
#    m2 = length(bo)
#    M = [zeros(Num, 3, 3) -AA'
#        AA zeros(Num, m1 + m2, m1 + m2)]
#    r = [qq; bb]
#    all_active_inds = collect(1:m1+m2)
#    Itr = powerset(all_active_inds) |> collect
#    Itr_reduced = []
#    for ass in Itr
#        if if_ass_feasible(ass, m1, m2)
#            push!(Itr_reduced, ass)
#        end
#    end
#    Itr = copy(Itr_reduced)

#    for active_inds in Itr
#        if !if_ass_feasible(active_inds, m1, m2)
#            continue
#        end
#        #length(active_inds) > 3 && continue
#        #assignment = [i ∈ active_inds for i in 1:m1+m2]
#        try
#            AA_active = collect(AA[active_inds, :])
#            bb_active = collect(bb[active_inds])
#            # TODO what if not unique primal? Need to resolve
#            if length(active_inds) == 3
#                zz = -AA_active \ bb_active
#                #zz = -(AA_active'*AA_active)\AA_active'*bb_active
#            else
#                #if linear system is underdetermined, use minimum norm solution (calculated by right inverse)
#                #Note: every solution has the same sd, i.e., zz[3], but different zz[1:2]
#                zz = -AA_active' * ((AA_active * AA_active') \ bb_active)
#            end
#            sd = zz[3]
#            sds[active_inds] = sd
#        catch err
#            if err isa LinearAlgebra.SingularException
#                continue
#            else
#                @warn(err)
#            end
#        end
#    end
#    sds
#end


function ConvexPolygon2DPointShrunk(P; c=0, s=0)
    A = P.A
    b = P.b
    V = P.V
    if c == 0
        c = sum(V) / length(V)
    end
    b = b + s * (A * c + b)
    PolyPlanning.ConvexPolygon2D(A, b)
end

function create_ass_playground(x0, ego_polys, obs_polys; fig=Figure(), ax=Axis(fig[1, 1], aspect=DataAspect()), θ=[], is_displaying=true)

    range_max = 10
    step_size = 0.1

    sg = SliderGrid(
        fig[2, 1],
        (label="x", range=-range_max:step_size:range_max, format="{:.1f}", startvalue=x0[1]),
        (label="y", range=-range_max:step_size:range_max, format="{:.1f}", startvalue=x0[2]),
        (label="θ", range=-range_max:step_size:range_max, format="{:.1f}", startvalue=x0[3]),
        (label="sd", range=-1.0:step_size:100.0, format="{:.1f}", startvalue=0),
    )

    for Pe in ego_polys
        for Po in obs_polys
            x = @lift([$(sg.sliders[1].value), $(sg.sliders[2].value), $(sg.sliders[3].value)])
            Aeb = @lift(PolyPlanning.shift_to(Pe.A, Pe.b, $x))
            Ae_shifted = GLMakie.lift(x -> x[1], Aeb)
            be_shifted = GLMakie.lift(x -> x[2], Aeb)
            Pe_shifted = @lift(PolyPlanning.ConvexPolygon2D($Ae_shifted, $be_shifted))

            # draw ego and obstacle
            PolyPlanning.plot!(ax, Pe_shifted; color=:blue, linewidth=3)
            PolyPlanning.plot!(ax, Po; color=:red, linewidth=3)

            # draw inflated ego and obstacle
            centroide = sum(Pe.V) / length(Pe.V)
            Ae = Pe.A
            be = Pe.b

            centroido = sum(Po.V) / length(Po.V)
            Ao = Po.A
            bo = Po.b
            primals_etc = @lift(get_single_sd_ids($x, Ae, be, centroide, Ao, bo, centroido))

            p_intercept_x = GLMakie.lift(x -> x[1][1], primals_etc)
            p_intercept_y = GLMakie.lift(x -> x[1][2], primals_etc)
            sd = GLMakie.lift(x -> x[1][3], primals_etc)

            ego_inflated = @lift(ConvexPolygon2DPointShrunk($Pe_shifted; c=0, s=$sd))
            obs_inflated = @lift(ConvexPolygon2DPointShrunk(Po; c=0, s=$sd))

            PolyPlanning.plot!(ax, ego_inflated; color=:lightblue, linewidth=2)
            PolyPlanning.plot!(ax, obs_inflated; color=:pink, linewidth=2)

            # draw arbitrary inflation
            sd_arbitrary = GLMakie.lift(x -> x, sg.sliders[4].value)

            ego_arb_inflated = @lift(ConvexPolygon2DPointShrunk($Pe_shifted; c=0, s=$sd_arbitrary))
            obs_arb_inflated = @lift(ConvexPolygon2DPointShrunk(Po; c=0, s=$sd_arbitrary))

            PolyPlanning.plot!(ax, ego_arb_inflated; color=:lightblue, linewidth=2, linestyle=:dash)
            PolyPlanning.plot!(ax, obs_arb_inflated; color=:pink, linewidth=2, linestyle=:dash)

            # draw the intercept and centroids
            R = @lift([cos($(x)[3]) sin($(x)[3])
                -sin($(x)[3]) cos($(x)[3])])
            centroidex = @lift($(x)[1:2] + $R * centroide)
            centroidex_x = GLMakie.lift(x -> x[1], centroidex)
            centroidex_y = GLMakie.lift(x -> x[2], centroidex)

            scatter!(ax, centroidex_x, centroidex_y; color=:blue)
            scatter!(ax, centroido[1], centroido[2]; color=:red)
            scatter!(ax, p_intercept_x, p_intercept_y; color=:green)

            # draw sds
            #sds_etc = @lift(g_col_single($x, Ae, be, centroide, Ao, bo, centroido))

            #AAbb = @lift(PolyPlanning.gen_LP_data($x, $Ae_shifted, $be_shifted, $centroidex, Ao, bo, centroido))
            #AA = GLMakie.lift(x -> x[1], AAbb)
            #bb = GLMakie.lift(x -> x[2], AAbb)

            ##@infiltrate
            #function filter_sds(sds_val, AA, bb)
            #    tol = 1e-4
            #    filtered = []
            #    for val in sds_val
            #        @infiltrate
            #        err = AA * val + bb
            #        if sum(err .>= -tol) == 8
            #            push!(filtered, val)
            #        end
            #    end
            #    filtered
            #end

            #sds_k = @lift(collect(keys($sds_etc)))
            #sds_val = @lift(collect(values($sds_etc)))
            #@infiltrate
            #@lift(filter_sds($sds_val, $AA, $bb))


            #point_text_obs = GLMakie.lift(sds_etc) do sds
            #    "Point: ($sds, $sds, $sds)"
            #end

            #text!(ax, point_text_obs, position=:lt, align=:left, fontsize=20, color=:black)

            #hss = map(1:8) do i
            #    HalfSpace(-AA[i, :], bb[i])
            #end
            #hr = hrep(hss) # H-representation for polyhedron 
            #poly = Polyhedra.polyhedron(hr)
            #vr = vrep(poly) # V-representation


            #mesh_poly = Polyhedra.Mesh(poly)
            #Makie.mesh(mesh_poly, color=:blue)
            #Makie.wireframe(mesh_poly)
        end
    end

    #scatter!(point, color=:red, markersize=20)

    limits!(ax, -range_max, range_max, -range_max, range_max)

    fig
end

Ve = [[-1, -1], [-1.1, 1.2], [0.6, 0.8], [1.7, -0.5]]
obs_polys = [PolyPlanning.ConvexPolygon2D(Ve)]

Vo = [[-1, -1], [-1, 1.0], [1, 1], [1, -1]]
ego_polys = [PolyPlanning.ConvexPolygon2D(Ve)]

x0 = [5.0, 0.0, 0.1, 0, 0, 0];
#obs_polys = PolyPlanning.gen_rect_obs(; a=0.5, b=2.0, x_shift=0.0);
#ego_polys = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);

create_ass_playground(x0, ego_polys, obs_polys)

#obs_polys = PolyPlanning.gen_rect_obs(; a=0.5, b=2.0, x_shift=0.0);
#ego_polys = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);
#
#Aes = [deepcopy(P.A) for P in ego_polys]
#bes = [deepcopy(P.b) for P in ego_polys]
#centroides = [sum(deepcopy(P.V)) / length(deepcopy(P.V)) for P in ego_polys]

#Ae = ego_polys[1].A
#be = ego_polys[1].b
#centroide = sum(ego_polys[1].V) / length(ego_polys[1].V)

#Ao = obs_polys[1].A
#bo = obs_polys[1].b
#centroido = sum(Po.V) / length(Po.V)

#assignments = get_single_sd_ids(xt, Ae, be, centroide, Ao, bo, centroido, derivs_per_sd, 0; is_newsd=is_newsd)

