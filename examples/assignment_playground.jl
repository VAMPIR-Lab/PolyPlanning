using PolyPlanning
using GLMakie
using JuMP
using Clp
using Combinatorics
using LinearAlgebra
using Polyhedra

# filter indices which are impossible to be active at the same time for one poly
function get_possible_constraint_ids(A, b)
    AA = Matrix(A)
    bb = b
    tol = 1e-3
    ind = collect(1:length(bb))
    inds = powerset(ind) |> collect
    itr = [i for i in inds if length(i)==2]
    feasible_inds=[]
    for i in itr
        try
            xx = - AA[i,:] \ bb[i]
            if all(AA*xx +bb .> -tol)
                push!(feasible_inds, i)
            end
        catch err
            if err isa LinearAlgebra.SingularException
                continue
            else
                # @warn(err)
            end
        end
    end
    feasible_inds
end

function get_single_sd_ids(xt, Ae, be, centroide, Ao, bo, centroido)
    Aex, bex = PolyPlanning.shift_to(Ae, be, xt)
    R = [cos(xt[3]) sin(xt[3])
        -sin(xt[3]) cos(xt[3])]
    centroidex = xt[1:2] + R * centroide
    AA, bb, qq = PolyPlanning.gen_LP_data(xt, Aex, bex, centroidex, Ao, bo, centroido)
    #@infiltrate
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

function if_ass_feasible(ass, m1, m2)
    if_ego = false
    if_obs = false
    for i in ass
        if i ∈ [i for i in 1:m1]
            if_ego = true
        end
        if i ∈ [i for i in m1+1:m1+m2]
            if_obs = true
        end
    end
    return length(ass) == 3 && if_ego && if_obs
end


function g_col_single(xt, Ae, be, centroide, Ao, bo, centroido; is_newsd=false)
    sds = Dict()
    Aex, bex = PolyPlanning.shift_to(Ae, be, xt)
    R = [cos(xt[3]) sin(xt[3])
        -sin(xt[3]) cos(xt[3])]
    centroidex = xt[1:2] + R * centroide
    AA, bb, qq = PolyPlanning.gen_LP_data(xt, Aex, bex, centroidex, Ao, bo, centroido)
    m1 = length(bex)
    m2 = length(bo)

    # get 32 assignments for 2 quadrilaterals
    inds_e = get_possible_constraint_ids(Ae, be)
    inds_o = get_possible_constraint_ids(Ao, bo)
    for i in 1:length(inds_o)
        inds_o[i] += [m1 , m1]
    end
    Itr=[]
    for i in 1:m1
        for ind in inds_o
            push!(Itr, sort(vcat(ind, i)))
        end
    end
    for i in m1+1:m1+m2
        for ind in inds_e
            push!(Itr, sort(vcat(ind, i)))
        end
    end

    # get 48 assignments for 2 quadrilaterals
    # all_active_inds = collect(1:m1+m2)
    # Itr = powerset(all_active_inds) |> collect
    # Itr_reduced = []
    # for ass in Itr
    #     if if_ass_feasible(ass, m1, m2)
    #         push!(Itr_reduced, ass)
    #     end
    # end
    # Itr = copy(Itr_reduced)
    for active_inds in Itr
        #length(active_inds) > 3 && continue
        #assignment = [i ∈ active_inds for i in 1:m1+m2]
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
    sds
end


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

function create_ass_playground(x0, ego_polys, obs_polys; fig=Figure(), θ=[], is_displaying=true)

    range_max = 10
    step_size = 0.1
    dim = size(ego_polys[1].A)[2]

    ax = Axis(fig[1, 2], aspect=DataAspect())
    # ax3 = Axis3(fig[2, 1])
    ax3 = LScene(fig[1, 1], scenekw = (camera = cam3d!, show_axis = true))
    sg = SliderGrid(
        fig[2, 1:2],
        (label="x", range=-range_max:step_size:range_max, format="{:.1f}", startvalue=x0[1]),
        (label="y", range=-range_max:step_size:range_max, format="{:.1f}", startvalue=x0[2]),
        (label="θ", range=-2π:π/100:2π, format="{:.2f}", startvalue=x0[3]),
        (label="sd", range=-1.0:step_size:100.0, format="{:.1f}", startvalue=0),
        (label="z lim", range=0.1:step_size:100.0, format="{:.1f}", startvalue=10.0)
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
            #@infiltrate
            sds_etc = @lift(g_col_single($x, Ae, be, centroide, Ao, bo, centroido))

            AAbb = @lift(PolyPlanning.gen_LP_data($x, $Ae_shifted, $be_shifted, $centroidex, Ao, bo, centroido))
            AA = GLMakie.lift(x -> x[1], AAbb)
            bb = GLMakie.lift(x -> x[2], AAbb)
            #@infiltrate

            function filter_sds(sds_etc, AA, bb)
                tol = 1e-4
                filtered_values = []
                filtered_keys = []
                for (key, val) in sds_etc
                    if all(AA * val + bb .>= -tol)
                        #@infiltrate
                        push!(filtered_values, val)
                        push!(filtered_keys, key)
                    end
                end
                filtered_values, filtered_keys
            end

            sds_k = @lift(collect(keys($sds_etc)))
            sds_val = @lift(collect(values($sds_etc)))
            filtered_ran = @lift(filter_sds($sds_etc, $AA, $bb))
            perm = @lift(sortperm($filtered_ran[1], by = x -> x[3]))
            filtered = @lift([($filtered_ran[1])[$perm], ($filtered_ran[2])[$perm]])
            
            max_intercepts = 12
            sigdigits = 2
            intercept_obs = Dict()

            # draw all intercepts
            map(1:max_intercepts) do i
                intercept_obs[i] = Observable([NaN; NaN; NaN])
                x = GLMakie.lift(x -> x[1], intercept_obs[i])
                y = GLMakie.lift(x -> x[2], intercept_obs[i])
                sd = GLMakie.lift(x -> x[3], intercept_obs[i])
                scatter!(ax, x, y; color=:yellow)
                text!(ax, x, y; align=(:center, :center), text=@lift("$(round($sd; sigdigits))"), color=:black, fontsize=10)
            end

            GLMakie.lift(filtered) do filtered
                #@infiltrate
                for (i, f) in enumerate(filtered[1])
                    if i > max_intercepts
                        @warn("$i too many intercepts")
                        continue
                    end
                    intercept_obs[i][] = copy(f[1:3])
                end

                for j in length(filtered[1])+1:max_intercepts
                    intercept_obs[j][] = [NaN, NaN, NaN]
                end
            end

            # print list
            point_text_obs = GLMakie.lift(filtered) do filtered
                if isempty(filtered[1])
                    return ""
                end
                sigdigits = 3
                #@infiltrate
                mapreduce(*, zip(filtered[1], filtered[2])) do (f, i)
                    "$(i), $(round(f[3]; sigdigits)), $(round.(f[1:2]; sigdigits) )\n"
                end
            end
            text!(ax, -5, -3; align=(:left, :top), text=point_text_obs, color=:black, fontsize=10)

            # manually add a bound to poly
            z_lim = GLMakie.lift(x -> x, sg.sliders[5].value)
            A_bound = zeros(dim+1)
            A_bound[end] = 1
            b_bound = @lift($primals_etc[1][3] + $z_lim)

            # 3d plot
            hss = GLMakie.lift(AA, bb) do AA, bb
                map(1:length(bb)) do i
                    HalfSpace(-AA[i, :], bb[i])
                end
            end
            hss = @lift(push!($hss, HalfSpace(A_bound, $b_bound)))
            hr = @lift(hrep($hss)) # H-representation for polyhedron 
            poly = @lift(Polyhedra.polyhedron($hr))
            mesh_poly = @lift(Polyhedra.Mesh($poly))
            # GLMakie.mesh!(ax3, mesh_poly, color=:blue, alpha=0.1)
            # GLMakie.wireframe!(ax3.scene, mesh_poly)
            GLMakie.mesh!(ax3.scene, mesh_poly, color=:green)

            # # plot 3d ego and obs
            hss_ego = GLMakie.lift(AA, bb) do AA, bb
                map(1:length(be)) do i
                    HalfSpace(-AA[i, :], bb[i])
                end
            end
            hss_ego = @lift(push!($hss_ego, HalfSpace(A_bound, $b_bound)))
            hr_ego = @lift(hrep($hss_ego)) # H-representation for polyhedron 
            poly_ego = @lift(Polyhedra.polyhedron($hr_ego))
            mesh_poly_ego = @lift(Polyhedra.Mesh($poly_ego))
            GLMakie.wireframe!(ax3.scene, mesh_poly_ego, color=:blue)

            hss_obs = GLMakie.lift(AA, bb) do AA, bb
                map(length(be)+1:length(bb)) do i
                    HalfSpace(-AA[i, :], bb[i])
                end
            end
            hss_obs = @lift(push!($hss_obs, HalfSpace(A_bound, $b_bound)))
            hr_obs = @lift(hrep($hss_obs)) # H-representation for polyhedron 
            poly_obs = @lift(Polyhedra.polyhedron($hr_obs))
            mesh_poly_obs = @lift(Polyhedra.Mesh($poly_obs))
            GLMakie.wireframe!(ax3.scene, mesh_poly_obs, color=:red)

            # draw all intercepts on 3d
            intercept3_obs = Dict()

            map(1:max_intercepts) do i
                intercept3_obs[i] = Observable([0.0; 0; -Inf])
                x = GLMakie.lift(x -> x[1], intercept3_obs[i])
                y = GLMakie.lift(x -> x[2], intercept3_obs[i])
                sd = GLMakie.lift(x -> x[3], intercept3_obs[i])
                scatter!(ax3, x, y, sd; color=:yellow)
                text!(ax3, x, y, sd; align=(:center, :center), text=@lift("$(round($sd; sigdigits))"), color=:black, fontsize=10)
            end

            GLMakie.lift(filtered) do filtered
                for (i, f) in enumerate(filtered[1])
                    intercept3_obs[i][] = copy(f[1:3])
                end
            end
        end
    end

    #scatter!(point, color=:red, markersize=20)

    # z_lim = GLMakie.lift(x -> x, sg.sliders[5].value)
    limits!(ax, -range_max, range_max, -range_max, range_max)
    # @lift(limits!(ax3, -range_max, range_max, -range_max, range_max, -1.0, $z_lim))

    fig
end

#Ve = [[-1, -1], [-1.1, 1.2], [0.6, 0.8], [1.7, -0.5]]
#obs_polys = [PolyPlanning.ConvexPolygon2D(Ve)]

#Vo = [[-1, -1], [-1, 1.0], [1, 1], [1, -1]]
#ego_polys = [PolyPlanning.ConvexPolygon2D(Ve)]

#x0 = [5.0, 0.0, 0.1, 0, 0, 0];
x0 = [2.0, 0.0, 2.89, 0, 0, 0];
obs_polys = PolyPlanning.gen_rect_obs(; a=0.5, b=2.0, x_shift=0.0);
ego_polys = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);

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

