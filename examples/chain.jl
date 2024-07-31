using PolyPlanning
using Polyhedra
using LinearAlgebra
using SparseArrays
using GLMakie
using Symbolics

Ve = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]
Vo = [[-1.0, -1, -1], [1, -1, -1], [0, 1, -1], [0, 0, 5]]
Pe = PolyPlanning.ConvexPolygon3D(Ve)
Po = PolyPlanning.ConvexPolygon3D(Vo)
ego_polys = [Pe]
obs_polys = [Po]

mrp = [1,2,3]/4
e, θ = PolyPlanning.axis_angle_from_mrp(mrp)
err = mrp - PolyPlanning.mrp_from_axis_angle(e, θ)
trans =zeros(3)+ [5,2,3]
x0 = [trans; mrp; zeros(6)]

T=20
n_x = 12
n_u = 6
n_xu = n_x + n_u
n_z = T * n_xu
n_ego = length(ego_polys)
n_obs = length(obs_polys)
# n_side_ego = length(ego_polys[1].b)
# n_side_obs = length(obs_polys[1].b)
n_side_ego = maximum([length(i.b) for i in ego_polys]) # just for the buffer
n_side_obs = maximum([length(i.b) for i in obs_polys]) # just for the buffer
n_dyn_cons = T * n_x
n_env_cons = T * n_xu
# combin_2_from_n = n::Int -> n * (n - 1) ÷ 2
# n_sds = (combin_2_from_n(n_side_ego) * n_side_obs + n_side_ego * combin_2_from_n(n_side_obs))

z = Symbolics.@variables(z[1:n_z])[1] |> Symbolics.scalarize
xt = Symbolics.@variables(xt[1:n_x])[1] |> Symbolics.scalarize

R = PolyPlanning.R_from_mrp(xt[4:6])
l = xt[1:3]

A_ego = collect(Pe.A)
b_ego = Pe.b
V_ego = Pe.V
centr_ego = Pe.c
m1 = length(b_ego)

A_obs = collect(Po.A)
b_obs = Po.b
V_obs = Po.V
centr_obs = Po.c

R0 = PolyPlanning.R_from_mrp(x0[4:6])
l0 = x0[1:3]
AA = [A_ego*R0' A_ego*centr_ego+b_ego
                A_obs    A_obs*centr_obs+b_obs]
bb = [b_ego-A_ego*R0'*l0
    b_obs]



i=1
j=1
k=1
get_A = Dict()
get_b = Dict()
ass = [2,3,4,6]




























function get_Num_0_matrix(n_row, n_col)
    return [Num(0) for _ in 1:n_row, _ in 1:n_col]
end

function get_aibi_wrt_xt(R, l, xt)
    # Symbolics of one row of original ego constraints
    Ai = Symbolics.@variables(Ai[1:3])[1] |> Symbolics.scalarize
    bi = Symbolics.@variables(bi)[1]

    # parameterized Ai and bi w.r.t. xt
    a1, a2, a3 = Ai' * R'
    b1 = bi - [a1, a2, a3]' *l

    # expression of gradient
    a1_wrt_xt = Symbolics.gradient(a1, xt; simplify=false)
    a2_wrt_xt = Symbolics.gradient(a2, xt; simplify=false)
    a3_wrt_xt = Symbolics.gradient(a3, xt; simplify=false)
    b1_wrt_xt = Symbolics.gradient(b1, xt; simplify=false)

    J_a1_wrt_xt = Symbolics.jacobian(a1_wrt_xt, xt; simplify=false)
    J_a2_wrt_xt = Symbolics.jacobian(a2_wrt_xt, xt; simplify=false)
    J_a3_wrt_xt = Symbolics.jacobian(a3_wrt_xt, xt; simplify=false)
    J_b1_wrt_xt = Symbolics.jacobian(b1_wrt_xt, xt; simplify=false)

    # function of gradient
    a1_wrt_xt_fun = Symbolics.build_function(a1_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    a2_wrt_xt_fun = Symbolics.build_function(a2_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    a3_wrt_xt_fun = Symbolics.build_function(a3_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    b1_wrt_xt_fun = Symbolics.build_function(b1_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    
    J_a1_wrt_xt_fun = Symbolics.build_function(J_a1_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    J_a2_wrt_xt_fun = Symbolics.build_function(J_a2_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    J_a3_wrt_xt_fun = Symbolics.build_function(J_a3_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    J_b1_wrt_xt_fun = Symbolics.build_function(J_b1_wrt_xt, xt, Ai, bi; expression=Val(false))[2]

    return [a1_wrt_xt_fun, a2_wrt_xt_fun, a3_wrt_xt_fun, b1_wrt_xt_fun, J_a1_wrt_xt_fun, J_a2_wrt_xt_fun, J_a3_wrt_xt_fun, J_b1_wrt_xt_fun]
end
aibi_wrt_xt_functions = get_aibi_wrt_xt(R, l, xt)





function get_Ab_ego_wrt_xt_fun(xt, A_ego, b_ego, aibi_wrt_xt_functions)
    a1_wrt_xt_fun, a2_wrt_xt_fun, a3_wrt_xt_fun, b1_wrt_xt_fun, J_a1_wrt_xt_fun, J_a2_wrt_xt_fun, J_a3_wrt_xt_fun, J_b1_wrt_xt_fun = aibi_wrt_xt_functions
    
    A_ego_wrt_xt = [Num[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for _ in 1:m1, _ in 1:3]
    b_ego_wrt_xt = [Num[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for _ in 1:m1]

    J_A_ego_wrt_xt = [get_Num_0_matrix(12, 12) for _ in 1:m1, _ in 1:3]
    J_b_ego_wrt_xt = [get_Num_0_matrix(12, 12) for _ in 1:m1]
    # replace symbolics Ai and bi with values of A_ego and b_ego, get expressions only including xt
    for k in 1:m1
        a1_wrt_xt_fun(A_ego_wrt_xt[k, 1], xt, A_ego[k,:], b_ego[k])
        a2_wrt_xt_fun(A_ego_wrt_xt[k, 2], xt, A_ego[k,:], b_ego[k])
        a3_wrt_xt_fun(A_ego_wrt_xt[k, 3], xt, A_ego[k,:], b_ego[k])
        b1_wrt_xt_fun(b_ego_wrt_xt[k], xt, A_ego[k,:], b_ego[k])

        J_a1_wrt_xt_fun(J_A_ego_wrt_xt[k, 1], xt, A_ego[k,:], b_ego[k])
        J_a2_wrt_xt_fun(J_A_ego_wrt_xt[k, 2], xt, A_ego[k,:], b_ego[k])
        J_a3_wrt_xt_fun(J_A_ego_wrt_xt[k, 3], xt, A_ego[k,:], b_ego[k])
        J_b1_wrt_xt_fun(J_b_ego_wrt_xt[k], xt, A_ego[k,:], b_ego[k])
    end

    get_A_ego_wrt_xt_fun = Symbolics.build_function(A_ego_wrt_xt, xt; expression=Val(false))[2]
    get_b_ego_wrt_xt_fun = Symbolics.build_function(b_ego_wrt_xt, xt; expression=Val(false))[2]
    get_J_A_ego_wrt_xt_fun = Symbolics.build_function(J_A_ego_wrt_xt, xt; expression=Val(false))[2]
    get_J_b_ego_wrt_xt_fun = Symbolics.build_function(J_b_ego_wrt_xt, xt; expression=Val(false))[2]
    return get_A_ego_wrt_xt_fun, get_b_ego_wrt_xt_fun, get_J_A_ego_wrt_xt_fun, get_J_b_ego_wrt_xt_fun
end

get_A_ego_wrt_xt_fun, get_b_ego_wrt_xt_fun, get_J_A_ego_wrt_xt_fun, get_J_b_ego_wrt_xt_fun = get_Ab_ego_wrt_xt_fun(xt, A_ego, b_ego, aibi_wrt_xt_functions)


# update A_wrt_xt_buffer and b_wrt_xt_buffer
function get_Ab_wrt_xt_buffer!(A_wrt_xt_buffer, b_wrt_xt_buffer, ass, A_ego_wrt_xt_buffer, b_ego_wrt_xt_buffer)
    for (k, ind) in enumerate(ass)
        if ind <= m1
            A_wrt_xt_buffer[k, 1:3] = @view A_ego_wrt_xt_buffer[ind, :]
            b_wrt_xt_buffer[k] = @view b_ego_wrt_xt_buffer[ind][1:end] # it returns a 0-dim view without [1:end]
        end
    end
end

# update J_A_wrt_xt_buffer and J_b_wrt_xt_buffer
function get_J_Ab_wrt_xt_buffer!(J_A_wrt_xt_buffer, J_b_wrt_xt_buffer, ass, J_A_ego_wrt_xt_buffer, J_b_ego_wrt_xt_buffer)
    for (k, ind) in enumerate(ass)
        if ind <= m1
            J_A_wrt_xt_buffer[k, 1:3] = @view J_A_ego_wrt_xt_buffer[ind, :]
            J_b_wrt_xt_buffer[k] = @view J_b_ego_wrt_xt_buffer[ind][1:end, 1:end] # it returns a 0-dim view without [1:end]
        end
    end
end




# do not use fill(), because every element points to the same vector zeros(12)
A_ego_wrt_xt_buffer = [zeros(12) for _ in 1:m1, _ in 1:3]
b_ego_wrt_xt_buffer = [zeros(12) for _ in 1:m1]
J_A_ego_wrt_xt_buffer = [zeros(12, 12) for _ in 1:m1, _ in 1:3]
J_b_ego_wrt_xt_buffer = [zeros(12, 12) for _ in 1:m1]

A_wrt_xt_buffer = [zeros(12) for _ in 1:4, _ in 1:4]
b_wrt_xt_buffer = [zeros(12) for _ in 1:4]
J_A_wrt_xt_buffer = [zeros(12, 12) for _ in 1:4, _ in 1:4]
J_b_wrt_xt_buffer = [zeros(12, 12) for _ in 1:4]


get_A_ego_wrt_xt_fun(A_ego_wrt_xt_buffer, x0)
get_b_ego_wrt_xt_fun(b_ego_wrt_xt_buffer, x0)
get_J_A_ego_wrt_xt_fun(J_A_ego_wrt_xt_buffer, x0)
get_J_b_ego_wrt_xt_fun(J_b_ego_wrt_xt_buffer, x0)


get_Ab_wrt_xt_buffer!(A_wrt_xt_buffer, b_wrt_xt_buffer, ass, A_ego_wrt_xt_buffer, b_ego_wrt_xt_buffer)
get_J_Ab_wrt_xt_buffer!(J_A_wrt_xt_buffer, J_b_wrt_xt_buffer, ass, J_A_ego_wrt_xt_buffer, J_b_ego_wrt_xt_buffer)



function get_sd_wrt_Ab_fun()
    # Symbolics of a 4d linear system
    A = Symbolics.@variables(A[1:4, 1:4])[1] |> Symbolics.scalarize
    b = Symbolics.@variables(b[1:4])[1] |> Symbolics.scalarize
    sd = (-A \ b)[4]

    # derivative of signed distance w.r.t. A and b
    sd_wrt_A = [Num(0) for _ in 1:4, _ in 1:4]
    sd_wrt_b = [Num(0) for _ in 1:4]
    J_sd_wrt_A = [Num(0) for _ in 1:4, _ in 1:4]
    J_sd_wrt_b = [Num(0) for _ in 1:4]
    for i in 1:4
        sd_wrt_b[i] = Symbolics.derivative(sd, b[i]; simplify=false)
        J_sd_wrt_b[i] = Symbolics.derivative(sd_wrt_b[i], b[i]; simplify=false)
        for j in 1:4
            sd_wrt_A[i, j] = Symbolics.derivative(sd, A[i, j]; simplify=false)
            J_sd_wrt_A[i, j] = Symbolics.derivative(sd_wrt_A[i, j], A[i, j]; simplify=false)
        end
    end

    get_sd_wrt_A_fun = Symbolics.build_function(sd_wrt_A, A, b; expression=Val(false))[2]
    get_sd_wrt_b_fun = Symbolics.build_function(sd_wrt_b, A, b; expression=Val(false))[2]
    get_J_sd_wrt_A_fun = Symbolics.build_function(J_sd_wrt_A, A, b; expression=Val(false))[2]
    get_J_sd_wrt_b_fun = Symbolics.build_function(J_sd_wrt_b, A, b; expression=Val(false))[2]
    return get_sd_wrt_A_fun, get_sd_wrt_b_fun, get_J_sd_wrt_A_fun, get_J_sd_wrt_b_fun
end
get_sd_wrt_A_fun, get_sd_wrt_b_fun, get_J_sd_wrt_A_fun, get_J_sd_wrt_b_fun = get_sd_wrt_Ab_fun()

sd_wrt_A_buffer = zeros(4,4)
sd_wrt_b_buffer = zeros(4)
J_sd_wrt_A_buffer = zeros(4,4)
J_sd_wrt_b_buffer = zeros(4)

get_sd_wrt_A_fun(sd_wrt_A_buffer, AA[ass,:], bb[ass])
get_sd_wrt_b_fun(sd_wrt_b_buffer, AA[ass,:], bb[ass])
get_J_sd_wrt_A_fun(J_sd_wrt_A_buffer, AA[ass,:], bb[ass])
get_J_sd_wrt_b_fun(J_sd_wrt_b_buffer, AA[ass,:], bb[ass])

# chain rule
sd_wrt_xt = sum(sd_wrt_A_buffer .* A_wrt_xt_buffer) + sum(sd_wrt_b_buffer .* b_wrt_xt_buffer)
J_sd_wrt_xt = sum(J_sd_wrt_A_buffer .* (A_wrt_xt_buffer * A_wrt_xt_buffer')) + sum(sd_wrt_A_buffer .* J_A_wrt_xt_buffer) + sum(J_sd_wrt_b_buffer .* (b_wrt_xt_buffer * b_wrt_xt_buffer')) + sum(sd_wrt_b_buffer .* J_b_wrt_xt_buffer)




# for i in 1:4
#     for j in 1:4
#         println(sum(der[i,j] - A_wrt_xt_buffer[i, j] * sd_wrt_A_buffer[i,j]))
#     end
# end

# test_J = [J_b_wrt_xt_buffer[1][1:6, 1:6], J_b_wrt_xt_buffer[2][1:6, 1:6], J_b_wrt_xt_buffer[3][1:6, 1:6], J_b_wrt_xt_buffer[4][1:6, 1:6]]

# t=time()
# for indi in 1:100000
#     sum(sd_wrt_b_buffer .* test_J)
# end
# println(time()-t)



# # test_J = [J_A_wrt_xt_buffer[x,y][4:6,4:6] for x in 1:4, y in 1:4]

# t=time()
# for indi in 1:100000
#     sd_wrt_b_buffer' * test_J
# end
# println(time()-t)

