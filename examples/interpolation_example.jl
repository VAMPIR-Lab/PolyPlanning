using Revise
using PolyPlanning
using GLMakie


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

# generate polys
# poly = PolyPlanning.gen_polys(2)
# P = poly[1]
# O = poly[2]
PV = [[-2.979186916878823, 3.195803611422573],
[-1.2947084695891702, 2.6586598643189405],
[0.0773098950420259, 0.30559870181761584],
[-0.3380371995350262, -1.8981621270724238]]
OV = [ [-3.1770398371954123, 3.3066601649153164],
[0.035741482233899724, 5.322223678322559],
[3.062821815219714, 1.5618522073600964],
[-1.27998284802841, 2.309513564931549]]
P = PolyPlanning.ConvexPolygon2D(PV)
O = PolyPlanning.ConvexPolygon2D(OV)
c = sum(P.V) / length(P.V)

iter = 1e6
eps = 1e-5
x1 = [-2., -2., 0.]
for i in 0:10
    delta = [10, 5, i* pi/10]
    x2 = x1 + delta
    @info i
    generate_s_t_plot(P, O, c, x1, x2; num_interp = 1000, interval_display = 100, eps_abs=eps, eps_rel=eps, max_iter=iter)#, max_iter=1e4)
end

# """
index = 314
x3 = x1 + (index-1)/999* delta
# x4 = x1 + .177* delta
P3 = PolyPlanning.shift_to_around_point(P, x3; c=c)
# P4 = PolyPlanning.shift_to_around_point(P, x4; c=c)
iter = 1e6
eps = 1e-5
plot_sd_point_shrunk(P3, O; eps_abs=eps, eps_rel=eps, max_iter=iter)
# plot_sd_point_shrunk(P4, O; eps_abs=1e-6, eps_rel=1e-6)
# """

LPdata = PolyPlanning.get_LP_data_point_shrunk(P3, O)
#PolyPlanning.get_single_sd_point_shrunk(P3, O;eps_abs = 1e-7, eps_rel=1e-7)
iter = [1e6, 1e7]
eps = [1e-4, 1e-5, 1e-6, 1e-7]
for i in iter
    for e in eps
        # @info "parameters" i=i e=e
        # plot_sd_point_shrunk(P3, O; eps_abs=e, eps_rel=e, max_iter=i)
    end
end

