using PolyPlanning
using GLMakie
using JuMP
using Clp
using Combinatorics
using LinearAlgebra
using Polyhedra

# filter indices which are impossible to be active at the same time for one poly
function get_3_possible_constraint_ids(A, b; tol=1e-4)
    AA = Matrix(A)
    dim = size(AA)[2]
    bb = b
    ind = collect(1:length(bb))
    inds = powerset(ind) |> collect
    itr = [i for i in inds if length(i)==dim]
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

# get possible pair of indices, which means edges
function get_2_possible_constraint_ids(A, b, V; tol=1e-4)
    matrix_v_face = map(V) do v
        A * v + b .< tol
    end
    feasible_inds = []
    for i in 1:length(matrix_v_face)
        for j in i+1:length(matrix_v_face)
            common_faces = matrix_v_face[i] + matrix_v_face[j] .== 2
            if sum(common_faces) == 2
                push!(feasible_inds, findall(common_faces))
            end
        end
    end
    feasible_inds
end

# get_possible_assignments_3d(Pe.A, Pe.b, Pe.V, Po.A, Po.b, Po.V)
# enumerate possible assignments (2 indices from one poly, and 2 indices from the other; or 3 + 1; or 1 + 3)
function get_possible_assignments_3d(Ae, be, Ve, Ao, bo, Vo)
    m1 = length(be)
    m2 = length(bo)
    inds_e = get_3_possible_constraint_ids(Ae, be)
    inds_o = get_3_possible_constraint_ids(Ao, bo)
    for i in eachindex(inds_o)
        inds_o[i] += [m1, m1, m1]
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

    inds_e = get_2_possible_constraint_ids(Ae, be, Ve)
    inds_o = get_2_possible_constraint_ids(Ao, bo, Vo)
    for i in eachindex(inds_o)
        inds_o[i] += [m1, m1]
    end
    for e in inds_e
        for o in inds_o
            push!(Itr, [e; o])
        end
    end

    Itr
end

function get_single_sd_ids(xt, Ae, be, centroide, Ao, bo, centroido)
    Aex, bex = PolyPlanning.shift_to(Ae, be, xt)
    R = [cos(xt[3]) sin(xt[3])
        -sin(xt[3]) cos(xt[3])]
    centroidex = xt[1:2] + R * centroide
    AA, bb, qq = PolyPlanning.gen_LP_data(Aex, bex, centroidex, Ao, bo, centroido)
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

function ConvexPolygon3DPointShrunk(P; c=0, s=0)
    A = P.A
    b = P.b
    V = P.V
    if c == 0
        c = sum(V) / length(V)
    end
    b = b + s * (A * c + b)
    PolyPlanning.ConvexPolygon3D(A, b)
end

function create_ass_playground(x0, ego_polys, obs_polys; fig=Figure(), θ=[], is_displaying=true)

    range_max = 10
    step_size = 0.1
    dim = size(ego_polys[1].A)[2]

    #ax = Axis(fig[1, 2], aspect=DataAspect())
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
            Pe_shifted = @lift(PolyPlanning.ConvexPolygon3D($Ae_shifted, $be_shifted))

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

            ego_inflated = @lift(ConvexPolygon3DPointShrunk($Pe_shifted; c=0, s=$sd))
            obs_inflated = @lift(ConvexPolygon3DPointShrunk(Po; c=0, s=$sd))

            PolyPlanning.plot_with_indices(ax, ego_inflated; color=:lightblue, linewidth=2)
            PolyPlanning.plot_with_indices(ax, obs_inflated; m1=length(be), color=:pink, linewidth=2)

            # draw arbitrary inflation
            sd_arbitrary = GLMakie.lift(x -> x, sg.sliders[4].value)

            ego_arb_inflated = @lift(ConvexPolygon3DPointShrunk($Pe_shifted; c=0, s=$sd_arbitrary))
            obs_arb_inflated = @lift(ConvexPolygon3DPointShrunk(Po; c=0, s=$sd_arbitrary))

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
            sds_etc = GLMakie.lift(x) do x
                sds, intercepts, AA, bb, sds_etc = PolyPlanning.g_col_single_3d(x, Ae, be, centroide, Ao, bo, centroido)
                sds_etc
            end

            AAbb = @lift(PolyPlanning.gen_LP_data($Ae_shifted, $be_shifted, $centroidex, Ao, bo, centroido))
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
            
            max_intercepts = 32
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

Ve = [[.25, -2, -1], [.25, 2, -1], [-.25, 2, -1], [-.25, -2, -1], [-.5, -.5, 1], [.5, -.5, 1], [-.5, .5, 1], [.5, .5, 1]]
Pe = PolyPlanning.ConvexPolygon3D(Ve)
ego_polys = [Pe]
Vo = [[-1.0, -1, -1], [1, -1, -1], [-1, 1, -1], [1, 1, -1], [0, 0, 5]]
Po = PolyPlanning.ConvexPolygon3D(Vo)
obs_polys = [Po]

mrp = [1,2,3]
e, θ = PolyPlanning.axis_angle_from_mrp(mrp)
err = mrp - PolyPlanning.mrp_from_axis_angle(e, θ)
if norm(err)>1e-4
    @warn err
end
println(e)
println(rad2deg(θ))
trans =zeros(3)+ [1,2,3]
x0 = [trans;mrp;zeros(6)]

Ae = Pe.A
be = Pe.b
Aet, bet = PolyPlanning.shift_to_3D(Ae, be, x0)
Pet = PolyPlanning.ConvexPolygon3D(Aet, bet)


fig = Figure()
ax3 = LScene(fig[1, 1], scenekw = (camera = cam3d!, show_axis = true))
PolyPlanning.plot_with_indices(ax3, Pet; color=:blue)
PolyPlanning.plot_with_indices(ax3, Po; m1=length(Pet.b), color=:red)
fig

get_possible_assignments_3d(Pe.A, Pe.b, Po.A, Po.b)

# fig = create_ass_playground(x0, ego_polys, obs_polys)
# GLMakie.save("./plots/playground.png", fig)
