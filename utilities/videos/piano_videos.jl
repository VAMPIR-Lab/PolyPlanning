using PolyPlanning
include("record_videos.jl")

ego_width = 0.5;
ego_length = 2.0;
ego_polys = PolyPlanning.gen_ego_rect(; a=ego_width, b=ego_length)
width = round(sqrt((ego_length / 2)^2 + ego_width^2); sigdigits=2) + 0.1
pre_L_length_base = 5.0
post_L_length = 3.0
init_x = pre_L_length_base - ego_length
init_y = -post_L_length - width / 2
pre_L_length = pre_L_length_base - width / 2
obs_polys = PolyPlanning.gen_L_corridor(; width, pre_L_length, post_L_length)
x0 = [init_x, init_y, 0.05, 0, 0, 0];
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
framerate=10

video_name = "piano"
record_videos(video_name, ego_polys, obs_polys, T, dt, R_cost, Q_cost, u1_max, u2_max, u3_max, xlims, ylims, framerate)
