function f(z, T, R)
    cost = 0.0
    for t in 1:T
        xt = @view(z[(t-1)*9+1:(t-1)*9+6])
        ut = @view(z[(t-1)*9+7:(t-1)*9+9])
        #cost += ut'*R*ut - goal_dir'*xt[1:2]
        #cost += xt[1:2]'*xt[1:2]
        #cost += -0.01*goal_dir'*xt[1:2] + ut'*R*ut
        #cost += 0.5*ut'*R*ut+ 0.001*xt'*xt
        cost += 0.5 * ut' * R * ut + 0.01 * xt[1:2]' * xt[1:2]
        #cost += 0.1*(xt[1:2]-goal_dir)'*(xt[1:2]-goal_dir) + ut'*R*ut
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
    x + dt * [x[4:6]; u[1:2]; u[3] / 10.0]
    #x + dt*u
end

function g_dyn(z, x0, T, dt, L)
    g = Num[]
    x_prev = x0
    for t in 1:T
        xt = @view(z[(t-1)*9+1:(t-1)*9+6])
        ut = @view(z[(t-1)*9+7:(t-1)*9+9])
        #append!(g, xt - kinematic_bicycle_dyn(x_prev, ut, dt, L))
        #append!(g, xt - pointmass_dyn(x_prev, ut, dt))
        append!(g, xt - identity_dyn(x_prev, ut, dt))
        x_prev = xt
    end
    g
end

function g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)
    g = Num[]
    for t in 1:T
        xt = @view(z[(t-1)*9+1:(t-1)*9+6])
        ut = @view(z[(t-1)*9+7:(t-1)*9+9])
        #append!(g, [p1_max + xt[1], p1_max-xt[1], xt[2]-p2_min, xt[3]-2π, 2π-xt[3]])
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

    #colors = [:red, :orange, :yellow, :green, :red, :orange]
    colors = [:red for _ in 1:length(polys)]
    for (P, c) in zip(polys, colors)
        plot!(ax, P; color=c)
    end
    display(fig)
end

function gen_gap(; width=1.25, length=0.25)
    l = length / 2
    w = width / 2
    os = l / 10 # offset
    xs = -3.0
    obs_len = 5

    p1 = PolyPlanning.ConvexPolygon2D([[-l + xs, w], [l + xs, w], [l + os + xs, obs_len], [-l - os + xs, obs_len]])
    p2 = PolyPlanning.ConvexPolygon2D([[-l + xs, -w], [l + xs, -w], [l + os + xs, -obs_len], [-l - os + xs, -obs_len]])

    [p1, p2]
end

function gen_polys(N; side_length=4)
    polys = map(1:N) do i
        offset = 1.5 * randn(2)
        P = PolyPlanning.ConvexPolygon2D([2 * randn(2) + offset for _ in 1:side_length])
        while length(P.V) != side_length
            P = PolyPlanning.ConvexPolygon2D([2 * randn(2) + offset for _ in 1:side_length])
        end
        P
    end
end

#  _______ a
# |______|  
# <-----> l*a
function gen_ego_rect(; a=0.5)
    offset = a / 20
    l_multip = 2
    p1 = ConvexPolygon2D([[0, 0], [0, a + offset], [l_multip * a - offset, a], [l_multip * a, 0]])
    [p1]
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
    p1 = ConvexPolygon2D([[0, 0], [0, a + offset], [l_multip * a - offset, a], [l_multip * a, 0]])
    p2 = ConvexPolygon2D([[0, a + offset], [0, h_multip * a], [a, h_multip * a], [a - offset, a]])
    [p1, p2]
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
