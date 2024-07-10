"""
Create a tuple of grouped index arrays. Assume the inner-most dimension changes most rapidly.

**usage:** 
```
for i in is
    for j in js 
        ...
        idxs[group_id][..., j, i]
    end
end
```

**e.g.:**
```
julia> idxs = get_sub2idxs((2,2),(3,4))
julia> idxs[1]
2×2 Matrix{Int64}:
 1  3
 2  4
julia> idxs[2]
3×4 Matrix{Int64}:
 5   8  11  14
 6   9  12  15
 7  10  13  16
```
"""
function get_sub2idxs(ns...)
    # reminder: julia uses col major matrices
    counter = 0
    idxs = map(ns) do n
        idx = zeros(Int, n...)
        prod_n = prod(n)
        idx[:] .= counter .+ collect(1:prod_n)
        counter += prod_n
        idx
    end
end

function f(z, T, R_cost, Q_cost)
    n_x = 6
    n_u = 3
    n_xu = n_x + n_u

    cost = 0.0
    for t in 1:T
        xt = @view(z[(t-1)*n_xu+1:(t-1)*n_xu+n_x])
        ut = @view(z[(t-1)*n_xu+n_x+1:(t-1)*n_xu+n_xu])
        cost += ut' * R_cost * ut + xt[1:2]' * Q_cost * xt[1:2]
    end
    cost
end

function pointmass_dyn(x, u, dt)
    p1, p2, v1, v2 = x
    a1, a2 = u
    x + dt * [v1 + dt / 2 * a1, v2 + dt / 2 * a2, a1, a2]
end

function kinematic_bicycle_dyn(x, u, dt, L)
    p1, p2, θ, v = x
    δ, a = u
    #x + dt * [cos(θ)*v, sin(θ)*v, tan(δ)*v / L, a]
    x + dt * [cos(θ + δ / 2) * (v + a / 2), sin(θ + δ / 2) * (v + a / 2), δ, a]
end

function identity_dyn(x, u, dt)
    #v = dt * [u[1:2]; u[3] / 10.0]
    #x + [dt * x[4:6] + dt * v / 2; v]
    x + dt * [x[4:6]; u[1:2]; u[3] / 10]
    #x + [u; u[1:2]; u[3] / 10.0]
    #x + dt*u
end

function g_dyn(z, x0, T, dt)
    n_x = 6
    n_u = 3
    n_xu = n_x + n_u

    g = Num[]
    x_prev = x0
    for t in 1:T
        xt = @view(z[(t-1)*n_xu+1:(t-1)*n_xu+n_x])
        ut = @view(z[(t-1)*n_xu+n_x+1:(t-1)*n_xu+n_xu])
        append!(g, xt - identity_dyn(x_prev, ut, dt))
        x_prev = xt
    end
    g
end

function g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)
    n_x = 6
    n_u = 3
    n_xu = n_x + n_u

    g = Num[]
    for t in 1:T
        xt = @view(z[(t-1)*n_xu+1:(t-1)*n_xu+n_x])
        ut = @view(z[(t-1)*n_xu+n_x+1:(t-1)*n_xu+n_xu])
        append!(g, [p1_max + xt[1], p1_max - xt[1], xt[2] - p2_min, u1_max - ut[1], ut[1] + u1_max, u2_max - ut[2], ut[2] + u2_max, u3_max - ut[3], ut[3] + u3_max,])
    end
    g
end

function shift_to(A, b, x)
    p = x[1:2]
    θ = x[3]
    R = [cos(θ) sin(θ)
        -sin(θ) cos(θ)]
    At = A * R'
    bt = b - At * p
    At, bt
end

function plot_polys(polys)
    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())

    colors = [:red for _ in 1:length(polys)]
    for (P, c) in zip(polys, colors)
        plot!(ax, P; color=c)
    end
    display(fig)
end

function gen_gap(; width=1.25, length=0.25, xs=-3.0)
    l = length / 2
    w = width / 2
    os = l / 10 # offset
    obs_len = 5.0

    p1 = PolyPlanning.ConvexPolygon2D([[-l + xs, w], [l + xs, w], [l + os + xs, obs_len], [-l - os + xs, obs_len]])
    p2 = PolyPlanning.ConvexPolygon2D([[-l + xs, -w], [l + xs, -w], [l + os + xs, -obs_len], [-l - os + xs, -obs_len]])

    [p1 p2]
end

## P1 
## ### P2
##
function gen_T_obstacle(; base_width=2.0, width=1.0, length=1.0)
    wall_width = 1.0
    #w = base_width / 2

    P1_top_left = [0, base_width / 2]
    P2_top_left = [wall_width, width / 2]

    P1 = ConvexPolygon2D([P1_top_left, P1_top_left + [wall_width, 0], P1_top_left + [wall_width, -base_width], P1_top_left + [0, -base_width]])
    P2 = ConvexPolygon2D([P2_top_left, P2_top_left + [length, 0], P2_top_left + [length, -width], P2_top_left + [0, -width]])

    [P1 P2]
end


#  ___    [  P1  ]
#| P2 |_______     | w
# ___| ___P3 |
function gen_L_corridor(; width=1.0, pre_L_length=1.0, post_L_length=1.0)
    wall_width = 4.0
    w = width / 2
    ext_multip = 2
    P1_top_left = [w, ext_multip * post_L_length]
    P2_top_left = [-wall_width - w, ext_multip * post_L_length]
    P3_top_left = [-w, -post_L_length - width]
    #P3_top_left_y =
    #a = 0.25
    offset = w / 20
    #l_multip = 5
    P1 = ConvexPolygon2D([P1_top_left, P1_top_left + [0, -(ext_multip + 1) * post_L_length], P1_top_left + [pre_L_length, -(ext_multip + 1) * post_L_length], P1_top_left + [pre_L_length + offset, offset]])
    P2 = ConvexPolygon2D([P2_top_left, P2_top_left + [-offset, -(ext_multip + 1) * post_L_length - width - wall_width - offset], P2_top_left + [wall_width, -(ext_multip + 1) * post_L_length - width - wall_width], P2_top_left + [wall_width, 0]])
    P3 = ConvexPolygon2D([P3_top_left, P3_top_left + [0, -wall_width], P3_top_left + [pre_L_length + width + offset, -wall_width - offset], P3_top_left + [pre_L_length + width, 0]])
    [P1 P2 P3]
end

# .| <- w -> | 
# ^origin
function gen_packing_wall(n_obs, n_sides; w=1.0, l=1.0, max_overlap=0.0)
    @assert n_sides > 2
    @assert n_obs > 0

    y_nom = 2 * l / n_obs

    polys = map(1:n_obs) do i
        y_min = -l + (i - 1) * y_nom
        y_max = -l + i * y_nom

        #base = [[0, y_min], [0, y_max]]
        base = [[0, y_min], [0, y_max + max_overlap * rand()]]

        prot = [[w * rand(), y_min - y_nom + (y_max - y_min + 2 * y_nom) * rand()] for _ in 1:n_sides-2]
        P = ConvexPolygon2D([base; prot])

        while length(P.V) != n_sides
            prot = [[w * rand(), y_min - y_nom + (y_max - y_min + 2 * y_nom) * rand()] for _ in 1:n_sides-2]
            P = ConvexPolygon2D([base; prot])
        end
        P
    end
end


function gen_polys(N; side_length=4)
    polys = map(1:N) do i
        offset = 3.5 * randn(2)
        P = PolyPlanning.ConvexPolygon2D([2 * randn(2) + offset for _ in 1:side_length])

        while length(P.V) != side_length
            P = PolyPlanning.ConvexPolygon2D([2 * randn(2) + offset for _ in 1:side_length])
        end

        function mean(v)
            sum(v) / length(v)
        end
        # move it closer to origin
        max_dist = 2
        dist_from_origin = sqrt(mean(P.V)'mean(P.V))
        while dist_from_origin >= max_dist
            Aeb = shift_to(P.A, P.b, [1.5 .* -mean(P.V) ./ dist_from_origin; 0])
            P = ConvexPolygon2D(Aeb[1], Aeb[2])
            dist_from_origin = sqrt(mean(P.V)'mean(P.V))
        end
        P
    end
end

function gen_rect_obs(; a=0.25, b=5 * a, x_shift=0.0)
    offset = 0 * a / 20
    P = ConvexPolygon2D([[x_shift, -b], [x_shift, b], [-a, b + offset], [-a, -b - offset]])
    [P]
end

function gen_simple_obs(; scale=1.0, offset=[0.0, 0.0])
    [ConvexPolygon2D([offset .+ [-1.0, -1], offset .+ [1.0, -1], offset .+ [1.0, 1], offset .+ [-1.0, 1]])]
end

#  _______ a
# |______|  
# <-----> b
function gen_ego_rect(; a=0.5, b=1.0)
    offset = 0 * a / 20
    P = ConvexPolygon2D([[0, 0], [0, a + offset], [b - offset, a], [b, 0]])

    # shift center of mass to the origin
    c = sum(P.V) / length(P.V)
    Ab = shift_to(P.A, P.b, [-c; 0])
    P = ConvexPolygon2D(Ab[1], Ab[2])
    [P]
end

# <-> a
# ---  ---
# | |__| |
# |______|  
# <-----> l*a
function gen_ego_U(; a=0.5)
    offset = a / 20
    l_multip = 5
    h_multip = 3
    p1 = ConvexPolygon2D([[0, 0], [0, a + offset], [l_multip * a - offset, a], [l_multip * a, 0]])
    p2 = ConvexPolygon2D([[0, a + offset], [0, h_multip * a], [a, h_multip * a], [a - offset, a]])
    p3 = ConvexPolygon2D([[(l_multip - 1) * a, a + offset], [(l_multip - 1) * a, h_multip * a], [l_multip * a, h_multip * a], [l_multip * a - offset, a]])
    [p1, p2, p3]
end
# x0 = [-5,.2,pi/2+.2,0,0,0]

# <-> a
# ---
# | |___
# |_____|  
# <-----> l*a
function gen_ego_L(; a=0.5)
    #p1 = ConvexPolygon2D([[-0.5, 0.55], [-0.55, 0], [1.0, 0], [1.0, 0.5]])
    #p2 = ConvexPolygon2D([[1.0, 0.5], [0.5, 0.45], [1.05, 1.5], [0.5, 1.5]])
    l_multip = 4
    h_multip = 4
    offset = a / 20
    P1 = ConvexPolygon2D([[0, 0], [0, a + offset], [l_multip * a - offset, a], [l_multip * a, 0]])
    P2 = ConvexPolygon2D([[0, a + offset], [0, h_multip * a], [a, h_multip * a], [a - offset, a]])

    # shift center of mass to the origin
    c = sum([P1.V; P2.V]) / length([P1.V; P2.V])
    Ab1 = shift_to(P1.A, P1.b, [-c; 0])
    Ab2 = shift_to(P2.A, P2.b, [-c; 0])
    P1 = ConvexPolygon2D(Ab1[1], Ab1[2])
    P2 = ConvexPolygon2D(Ab2[1], Ab2[2])

    [P1, P2]
end
# x0 = [-5.5,0,pi/2,0,0,0]

function gen_hallway()
    P1 = PolyPlanning.ConvexPolygon2D([[-5.0, 0], [-5, -2], [-0.4, -2], [-0.4, 0]])
    #P1 = PolyPlanning.ConvexPolygon2D([[-5.0,0],[-5,-1], [-1.4,-1], [-1.4,0]]);  #uncomment to make hallwawayy bigger
    P2 = PolyPlanning.ConvexPolygon2D([[5.0, 0], [5, -2.8], [0.4, -2.8], [0.4, 0]])
    P3 = PolyPlanning.ConvexPolygon2D([[-5.0, -2.8], [5, -2.8], [-5, -8], [5, -8]])
    P4 = PolyPlanning.ConvexPolygon2D([[-5.0, 0], [-5, -8], [-6, -8], [-6, 0]])
    [P1, P2, P3, P4]
end
