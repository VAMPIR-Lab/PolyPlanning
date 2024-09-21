using PolyPlanning
include("record_videos.jl")

x0 = [4., 0, pi / 2, 0, 0, 0];
ego_polys = PolyPlanning.gen_ego_L();
obs_polys = PolyPlanning.gen_gap(; width=1.5, xs=2.)

T = 20;
dt = 0.2;
R_cost = 1e-4 * PolyPlanning.I(3);
R_cost[3, 3] = R_cost[3, 3] / 100.0;
Q_cost = 5e-3 * PolyPlanning.I(2);
u1_max = 10.0;
u2_max = 10.0;
u3_max = π;

xlims = [-2, 6.]
ylims = [-4., 4.]
framerate=10

video_name = "L_piano"
record_videos(video_name, ego_polys, obs_polys, T, dt, R_cost, Q_cost, u1_max, u2_max, u3_max, xlims, ylims, framerate)
