#using PolyPlanning

#x0 = [-5.5, 0, pi / 2, 0, 0, 0];
#ego_polys = PolyPlanning.gen_ego_L();
#gap_polys = PolyPlanning.gen_gap();
#R_cost = 1e-3 * PolyPlanning.I(3);
#R_cost[3, 3] = R_cost[3, 3] / 100.0;
#T = 20;

#nonsmooth_prob = PolyPlanning.setup_nonsmooth(
#    ego_polys,
#    gap_polys;
#    T,
#    dt=0.2,
#    R_cost,
#    Q_cost=2e-3 * PolyPlanning.I(2),
#    u1_max=10.0,
#    u2_max=10.0,
#    u3_max=π,
#    n_sd_slots=2
#);

#our_sol = PolyPlanning.solve_nonsmooth(nonsmooth_prob, x0; is_displaying=true, sleep_duration=0.3, is_recording=true)


#(fig, update_fig, ax) = PolyPlanning.visualize_nonsmooth(x0, T, ego_polys, gap_polys)
#framerate=10
#hist_len = length(our_sol.θ_history)
#frames = 2:hist_len+1*framerate

#PolyPlanning.record(fig, "anim.mp4", frames;
#    framerate) do frame
#    #@infiltrate
#    if frame <= hist_len
#        update_fig(our_sol.θ_history[frame])
#    end
#end

using PolyPlanning

x0 = [5.0, 0.0, 0.1, 0, 0, 0];
obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);
ego_rect = PolyPlanning.gen_ego_rect(; a=0.5, b=1.0);
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;

sep_prob = PolyPlanning.setup_sep_planes(
    ego_rect,
    obs_polys;
    T=20,
    dt=0.2,
    Rf,
    Qf=5e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π
)

sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0; is_displaying=false, sleep_duration=0.25, is_recording=true)

(fig, update_fig, ax) = PolyPlanning.visualize_sep_planes(x0, T, ego_rect, obs_polys)

PolyPlanning.xlims!(ax, -2., 2.)
PolyPlanning.ylims!(ax, -2., 2.)
framerate=10
hist_len = length(sep_sol.θ_history)
frames = 2:hist_len+1*framerate

PolyPlanning.record(fig, "anim.mp4", frames;
    framerate) do frame
    #@infiltrate
    if frame <= hist_len
        update_fig(sep_sol.θ_history[frame])
    end
end