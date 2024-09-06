using PolyPlanning
using JLD2
using Dates
using GLMakie
using LinearAlgebra
using ProgressMeter
using PATHSolver
# experiment parameters (ignored if is_loading_exp or is_loading_res)
n_maps = 10 # number of maps
n_x0s = 100 # number of initial conditions
r=rand

Ve = [[-1-r(), -1-r(), -1-r()], [1+r(), -1-r(), -1-r()], [0, 1+r(), -1-r()], [0, 0, 1+r()]]
Pe = PolyPlanning.ConvexPolygon3D(Ve)
ego_polys = [Pe]

x0s = map(1:n_x0s) do i
    axis = rand(3)
    mrp = axis / norm(axis) * (0.8+0.4*r()) # norm(mrp)∈[0.8, 1.2]
    trans =zeros(3) + [5,2,3] + r(3)
    [trans; mrp]
end

maps = map(1:n_maps) do i
    Vo = [[-1-r(), -1-r(), -1-r()], [1+r(), -1-r(), -1-r()], [0, 1+r(), -1-r()], [0, 0, 5r()]]
    Po = PolyPlanning.ConvexPolygon3D(Vo)
    [Po]
end

our_sols = Dict()

p = Progress(n_maps * n_x0s, dt=1.0)
for (imap, obs_polys) in enumerate(maps)
    nonsmooth_prob = PolyPlanning.setup_nonsmooth_3d_0step(
        ego_polys,
        obs_polys;
        T=1,
        Q_cost=2e-3 * PolyPlanning.I(3),
        n_sd_slots=15
    )
    for (ix, x0) in enumerate(x0s)
        our_sol = PolyPlanning.solve_nonsmooth_3d_0step(nonsmooth_prob, x0; is_displaying=false)
        our_sols[imap, ix] = our_sol
        next!(p)
    end
end

success_rates_our = []
for imap in 1:n_maps
    map = fill(false, n_x0s)
    for ix in 1:n_x0s
        map[ix] = our_sols[imap, ix].status == PATHSolver.MCP_Solved
    end
    push!(success_rates_our, map)
end
@info success_rate = sum(sum(success_rates_our))/(n_maps*n_x0s)

