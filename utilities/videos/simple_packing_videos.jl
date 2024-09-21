using PolyPlanning

x0 = [2.0, 0, pi / 4, 0, 0, 0];
ego_width = 0.5;
ego_length = 2.0;
ego_polys = PolyPlanning.gen_ego_rect(; a=ego_width, b=ego_length);
obs_polys = PolyPlanning.gen_rect_obs(; a=0.25);
R_cost = 1e-3 * PolyPlanning.I(3);
R_cost[3, 3] = R_cost[3, 3] / 100.0;
T = 20;
dt = 0.2;
R_cost = 1e-3 * PolyPlanning.I(3);
R_cost[3, 3] = R_cost[3, 3] / 100.0;
Q_cost = 2e-3 * PolyPlanning.I(2);
u1_max = 10.0;
u2_max = 10.0;
u3_max = π;

xlims = [-1., 3.]
ylims = [-1.5, 1.5]
framerate=10

video_name = "simple_packing"
include("record_videos.jl")