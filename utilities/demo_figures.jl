using Dates
using PolyPlanning
using GLMakie
using CairoMakie 

exp_file_date = "2024-07-04_0428"
res_file_date = "2024-07-04_0428"
data_dir = "data"

plot_n_rows = 5
plot_n_cols = 5

# simple packing
exp_name = "simple_packing"
ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.load_all(exp_name, exp_file_date, res_file_date; data_dir)

GLMakie.activate!()
#PolyPlanning.visualize_multi(x0s, maps, our_sols, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
x0_idx = 16
map_idx = 1
x0 = x0s[x0_idx]
map = maps[map_idx]
sol = our_sols[(map_idx, x0_idx)]

# plot
fig=PolyPlanning.Figure()
ax=PolyPlanning.Axis(fig[1, 1], aspect=PolyPlanning.DataAspect())
PolyPlanning.hidedecorations!(ax)

CairoMakie.activate!()
(fig, update_fig, ax) = PolyPlanning.visualize_nonsmooth(x0, param.T, ego_poly, map; fig, ax,sol.res.θ, is_displaying=false)
GLMakie.save("./figures/simple_packing.pdf", fig) 

# simple gap
exp_name = "simple_gap"
ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.load_all(exp_name, exp_file_date, res_file_date; data_dir)

GLMakie.activate!()
#PolyPlanning.visualize_multi(x0s, maps, our_sols, param.T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
x0_idx = 106 # 56, 106, 186
map_idx = 2 # 5, 2,  3
x0 = x0s[x0_idx]
map = maps[map_idx]
sol = our_sols[(map_idx, x0_idx)]

# plot
fig=PolyPlanning.Figure()
ax=PolyPlanning.Axis(fig[1, 1], aspect=PolyPlanning.DataAspect())
PolyPlanning.hidedecorations!(ax)

CairoMakie.activate!()
(fig, update_fig, ax) = PolyPlanning.visualize_nonsmooth(x0, param.T, ego_poly, map; fig, ax,sol.res.θ, is_displaying=false)
GLMakie.save("./figures/simple_gap.pdf", fig) 

# piano
exp_name = "piano"
ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.load_all(exp_name, exp_file_date, res_file_date; data_dir)

GLMakie.activate!()
#PolyPlanning.visualize_multi(x0s, maps, our_sols, param.T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
x0_idx = 3 
map_idx = 1 
x0 = x0s[x0_idx]
map = maps[map_idx]
sol = our_sols[(map_idx, x0_idx)]

# plot
fig=PolyPlanning.Figure()
ax=PolyPlanning.Axis(fig[1, 1], aspect=PolyPlanning.DataAspect())
PolyPlanning.hidedecorations!(ax)

CairoMakie.activate!()
(fig, update_fig, ax) = PolyPlanning.visualize_nonsmooth(x0, param.T, ego_poly, map; fig, ax,sol.res.θ, is_displaying=false)
GLMakie.save("./figures/piano.pdf", fig) 

# random packing
exp_name = "random_packing"
ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.load_all(exp_name, exp_file_date, res_file_date; data_dir)

GLMakie.activate!()
#PolyPlanning.visualize_multi(x0s, maps, our_sols, param.T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
x0_idx = 13  # 11, 13, 44, 66, 71
map_idx = 5 # 5, 5, 5
x0 = x0s[x0_idx]
map = maps[map_idx]
sol = our_sols[(map_idx, x0_idx)]

# plot
fig=PolyPlanning.Figure()
ax=PolyPlanning.Axis(fig[1, 1], aspect=PolyPlanning.DataAspect())
PolyPlanning.hidedecorations!(ax)

CairoMakie.activate!()
(fig, update_fig, ax) = PolyPlanning.visualize_nonsmooth(x0, param.T, ego_poly, map; fig, ax,sol.res.θ, is_displaying=false)
GLMakie.save("./figures/random_packing.pdf", fig) 

# L through gap
exp_name = "L_piano"
ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.load_all(exp_name, exp_file_date, res_file_date; data_dir)

GLMakie.activate!()
#PolyPlanning.visualize_multi(x0s, maps, our_sols, param.T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
x0_idx = 63 # 56, 63
map_idx = 4 # 4, 4
x0 = x0s[x0_idx]
map = maps[map_idx]
sol = our_sols[(map_idx, x0_idx)]

# plot
fig=PolyPlanning.Figure()
ax=PolyPlanning.Axis(fig[1, 1], aspect=PolyPlanning.DataAspect())
PolyPlanning.hidedecorations!(ax)

CairoMakie.activate!()
(fig, update_fig, ax) = PolyPlanning.visualize_nonsmooth(x0, param.T, ego_poly, map; fig, ax,sol.res.θ, is_displaying=false)
GLMakie.save("./figures/L_piano.pdf", fig) 

# random L packing
exp_name = "random_L_packing"
ego_poly, x0s, maps, param = PolyPlanning.load_experiment(exp_name, exp_file_date; data_dir)
our_sols, sep_sols, dcol_sols, kkt_sols = PolyPlanning.load_all(exp_name, exp_file_date, res_file_date; data_dir)

GLMakie.activate!()
#PolyPlanning.visualize_multi(x0s, maps, our_sols, param.T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
x0_idx = 61 # 22, 61, 68
map_idx = 4 # 4,
x0 = x0s[x0_idx]
map = maps[map_idx]
sol = our_sols[(map_idx, x0_idx)]

# plot
fig=PolyPlanning.Figure()
ax=PolyPlanning.Axis(fig[1, 1], aspect=PolyPlanning.DataAspect())
PolyPlanning.hidedecorations!(ax)

CairoMakie.activate!()
(fig, update_fig, ax) = PolyPlanning.visualize_nonsmooth(x0, param.T, ego_poly, map; fig, ax,sol.res.θ, is_displaying=false)
GLMakie.save("./figures/random_L_packing.pdf", fig) 

GLMakie.activate!()