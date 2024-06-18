using PolyPlanning

x0 = [5.0, 0.0, 0.1, 0, 0, 0];
obs_polys = PolyPlanning.gen_rect_obs(; a=0.5, b=2.0);
ego_polys = PolyPlanning.gen_ego_rect(; a=0.5, b=2.0);

PolyPlanning.setup_nonsmooth(ego_polys, obs_polys)