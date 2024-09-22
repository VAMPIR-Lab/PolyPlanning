using PolyPlanning
include("record_videos.jl")

x0 = [5.0, 0.0, -π / 4, 0, 0, 0];
obs_polys = PolyPlanning.gen_packing_wall(4, 4; w=5.0, l=5.0, max_overlap=0.0, seed=42);
#PolyPlanning.plot_polys(obs_polys);
ego_polys = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);
T = 20;
dt = 0.2;
R_cost = 1e-4 * PolyPlanning.I(3);
R_cost[3, 3] = R_cost[3, 3] / 100.0;
Q_cost = 5e-3 * PolyPlanning.I(2);
u1_max = 10.0;
u2_max = 10.0;
u3_max = π;

xlims = [-.5, 6.5]
ylims = [-3.5, 3.5]
framerate=90

video_name = "random_packing"
record_videos(video_name, ego_polys, obs_polys, T, dt, R_cost, Q_cost, u1_max, u2_max, u3_max, xlims, ylims, framerate)
