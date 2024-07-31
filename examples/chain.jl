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
λsd = Symbolics.@variables(λsd)[1]

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
ass = [1,4,5,8]
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

function get_Ab_ego_wrt_xt_fun(xt, A_ego, b_ego, aibi_wrt_xt_functions)
    a1_wrt_xt_fun, a2_wrt_xt_fun, a3_wrt_xt_fun, b1_wrt_xt_fun, J_a1_wrt_xt_fun, J_a2_wrt_xt_fun, J_a3_wrt_xt_fun, J_b1_wrt_xt_fun = aibi_wrt_xt_functions
    
    A_ego_wrt_xt = [Num[0, 0, 0, 0, 0, 0] for _ in 1:m1, _ in 1:3]
    b_ego_wrt_xt = [Num[0, 0, 0, 0, 0, 0] for _ in 1:m1]

    J_A_ego_wrt_xt = [get_Num_0_matrix(6, 6) for _ in 1:m1, _ in 1:3]
    J_b_ego_wrt_xt = [get_Num_0_matrix(6, 6) for _ in 1:m1]
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


# update A_wrt_xt and b_wrt_xt
function get_Ab_wrt_xt!(A_wrt_xt, b_wrt_xt, ass, A_ego_wrt_xt, b_ego_wrt_xt)
    for (k, ind) in enumerate(ass)
        if ind <= m1
            A_wrt_xt[k, 1:3] .= A_ego_wrt_xt[ind, :]
            b_wrt_xt[k] .= b_ego_wrt_xt[ind]
        end
    end
end

# update J_A_wrt_xt and J_b_wrt_xt
function get_J_Ab_wrt_xt!(J_A_wrt_xt, J_b_wrt_xt, ass, J_A_ego_wrt_xt, J_b_ego_wrt_xt)
    for (k, ind) in enumerate(ass)
        if ind <= m1
            J_A_wrt_xt[k, 1:3] .= J_A_ego_wrt_xt[ind, :]
            J_b_wrt_xt[k] .= J_b_ego_wrt_xt[ind]
        end
    end
end

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
    # the last row is always from obs, the last column is also stable
    for i in 1:3
        sd_wrt_b[i] = Symbolics.derivative(sd, b[i]; simplify=false)
        J_sd_wrt_b[i] = Symbolics.derivative(sd_wrt_b[i], b[i]; simplify=false)
        for j in 1:3
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

function get_sd_wrt_xt!(sd_wrt_xt, sd_wrt_A, A_wrt_xt, sd_wrt_b, b_wrt_xt)
    # sd_wrt_xt .= 0
    # sd_wrt_xt .+= sum(sd_wrt_A .* A_wrt_xt) + sum(sd_wrt_b .* b_wrt_xt)
    sd_wrt_xt .= sum(sd_wrt_A .* A_wrt_xt) + sum(sd_wrt_b .* b_wrt_xt)
end

function get_J_sd_wrt_xt!(J_sd_wrt_xt, sd_wrt_A, A_wrt_xt, sd_wrt_b, b_wrt_xt, J_sd_wrt_A, J_A_wrt_xt, J_sd_wrt_b, J_b_wrt_xt)
    # J_sd_wrt_xt .= 0
    # J_sd_wrt_xt .+= sum(J_sd_wrt_A .* (A_wrt_xt * A_wrt_xt')) + sum(sd_wrt_A .* J_A_wrt_xt) + sum(J_sd_wrt_b .* (b_wrt_xt * b_wrt_xt')) + sum(sd_wrt_b .* J_b_wrt_xt)
    J_sd_wrt_xt .= sum(J_sd_wrt_A .* (A_wrt_xt * A_wrt_xt')) + sum(sd_wrt_A .* J_A_wrt_xt) + sum(J_sd_wrt_b .* (b_wrt_xt * b_wrt_xt')) + sum(sd_wrt_b .* J_b_wrt_xt)
end










get_A_ego_wrt_xt_fun = Dict()
get_b_ego_wrt_xt_fun = Dict()
get_J_A_ego_wrt_xt_fun = Dict()
get_J_b_ego_wrt_xt_fun = Dict()


# xt[7:12] don't show up
aibi_wrt_xt_functions = get_aibi_wrt_xt(R, l, xt[1:6])
get_sd_wrt_A_fun, get_sd_wrt_b_fun, get_J_sd_wrt_A_fun, get_J_sd_wrt_b_fun = get_sd_wrt_Ab_fun()
get_A_ego_wrt_xt_fun[i], get_b_ego_wrt_xt_fun[i], get_J_A_ego_wrt_xt_fun[i], get_J_b_ego_wrt_xt_fun[i] = get_Ab_ego_wrt_xt_fun(xt[1:6], A_ego, b_ego, aibi_wrt_xt_functions)












# do not use fill(), because every element points to the same vector zeros(6)
A_ego_wrt_xt_buffer = [zeros(6) for _ in 1:m1, _ in 1:3]
b_ego_wrt_xt_buffer = [zeros(6) for _ in 1:m1]
J_A_ego_wrt_xt_buffer = [zeros(6, 6) for _ in 1:m1, _ in 1:3]
J_b_ego_wrt_xt_buffer = [zeros(6, 6) for _ in 1:m1]

A_wrt_xt_buffer = [zeros(6) for _ in 1:4, _ in 1:4]
b_wrt_xt_buffer = [zeros(6) for _ in 1:4]
J_A_wrt_xt_buffer = [zeros(6, 6) for _ in 1:4, _ in 1:4]
J_b_wrt_xt_buffer = [zeros(6, 6) for _ in 1:4]

sd_wrt_A_buffer = zeros(4,4)
sd_wrt_b_buffer = zeros(4)
J_sd_wrt_A_buffer = zeros(4,4)
J_sd_wrt_b_buffer = zeros(4)

sd_wrt_xt_buffer = zeros(6)
J_sd_wrt_xt_buffer = zeros(6, 6)


sd_lag_buf = zeros(n_x)
Jsdlag_buf = zeros(n_x, n_x+1)
Jsd_buf = zeros(n_x)





# get_A_ego_wrt_xt_fun[i](A_ego_wrt_xt_buffer, x0[1:6])
# get_b_ego_wrt_xt_fun[i](b_ego_wrt_xt_buffer, x0[1:6])
# get_J_A_ego_wrt_xt_fun[i](J_A_ego_wrt_xt_buffer, x0[1:6])
# get_J_b_ego_wrt_xt_fun[i](J_b_ego_wrt_xt_buffer, x0[1:6])


# get_Ab_wrt_xt!(A_wrt_xt_buffer, b_wrt_xt_buffer, ass, A_ego_wrt_xt_buffer, b_ego_wrt_xt_buffer)
# get_J_Ab_wrt_xt!(J_A_wrt_xt_buffer, J_b_wrt_xt_buffer, ass, J_A_ego_wrt_xt_buffer, J_b_ego_wrt_xt_buffer)


# get_sd_wrt_A_fun(sd_wrt_A_buffer, AA[ass,:], bb[ass])
# get_sd_wrt_b_fun(sd_wrt_b_buffer, AA[ass,:], bb[ass])
# get_J_sd_wrt_A_fun(J_sd_wrt_A_buffer, AA[ass,:], bb[ass])
# get_J_sd_wrt_b_fun(J_sd_wrt_b_buffer, AA[ass,:], bb[ass])


# # chain rule
# sd_wrt_xt = sum(sd_wrt_A_buffer .* A_wrt_xt_buffer) + sum(sd_wrt_b_buffer .* b_wrt_xt_buffer)
# J_sd_wrt_xt = sum(J_sd_wrt_A_buffer .* (A_wrt_xt_buffer * A_wrt_xt_buffer')) + sum(sd_wrt_A_buffer .* J_A_wrt_xt_buffer) + sum(J_sd_wrt_b_buffer .* (b_wrt_xt_buffer * b_wrt_xt_buffer')) + sum(sd_wrt_b_buffer .* J_b_wrt_xt_buffer)

# get_sd_wrt_xt!(sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer)
# get_J_sd_wrt_xt!(J_sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer, J_sd_wrt_A_buffer, J_A_wrt_xt_buffer, J_sd_wrt_b_buffer, J_b_wrt_xt_buffer)







# Jsd = ∇sd, 12d vector
function get_Jsd!(J_sd, xt, i, ass, AA, bb)
    get_A_ego_wrt_xt_fun[i](A_ego_wrt_xt_buffer, xt[1:6])
    get_b_ego_wrt_xt_fun[i](b_ego_wrt_xt_buffer, xt[1:6])
    # get_J_A_ego_wrt_xt_fun[i](J_A_ego_wrt_xt_buffer, xt[1:6])
    # get_J_b_ego_wrt_xt_fun[i](J_b_ego_wrt_xt_buffer, xt[1:6])

    get_Ab_wrt_xt!(A_wrt_xt_buffer, b_wrt_xt_buffer, ass, A_ego_wrt_xt_buffer, b_ego_wrt_xt_buffer)
    # get_J_Ab_wrt_xt!(J_A_wrt_xt_buffer, J_b_wrt_xt_buffer, ass, J_A_ego_wrt_xt_buffer, J_b_ego_wrt_xt_buffer)

    get_sd_wrt_A_fun(sd_wrt_A_buffer, AA[ass,:], bb[ass])
    get_sd_wrt_b_fun(sd_wrt_b_buffer, AA[ass,:], bb[ass])
    # get_J_sd_wrt_A_fun(J_sd_wrt_A_buffer, AA[ass,:], bb[ass])
    # get_J_sd_wrt_b_fun(J_sd_wrt_b_buffer, AA[ass,:], bb[ass])

    get_sd_wrt_xt!(sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer)
    # get_J_sd_wrt_xt!(J_sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer, J_sd_wrt_A_buffer, J_A_wrt_xt_buffer, J_sd_wrt_b_buffer, J_b_wrt_xt_buffer)
    
    J_sd .= [sd_wrt_xt_buffer; zeros(6)]
end

get_Jsd!(Jsd_buf, x0, i, ass, AA, bb)

# sd_lag = -λsd*∇sd, 12d vector
function get_sd_lag!(sd_lag, xt, λsd, i, ass, AA, bb)
    get_A_ego_wrt_xt_fun[i](A_ego_wrt_xt_buffer, xt[1:6])
    get_b_ego_wrt_xt_fun[i](b_ego_wrt_xt_buffer, xt[1:6])
    # get_J_A_ego_wrt_xt_fun[i](J_A_ego_wrt_xt_buffer, xt[1:6])
    # get_J_b_ego_wrt_xt_fun[i](J_b_ego_wrt_xt_buffer, xt[1:6])

    get_Ab_wrt_xt!(A_wrt_xt_buffer, b_wrt_xt_buffer, ass, A_ego_wrt_xt_buffer, b_ego_wrt_xt_buffer)
    # get_J_Ab_wrt_xt!(J_A_wrt_xt_buffer, J_b_wrt_xt_buffer, ass, J_A_ego_wrt_xt_buffer, J_b_ego_wrt_xt_buffer)

    get_sd_wrt_A_fun(sd_wrt_A_buffer, AA[ass,:], bb[ass])
    get_sd_wrt_b_fun(sd_wrt_b_buffer, AA[ass,:], bb[ass])
    # get_J_sd_wrt_A_fun(J_sd_wrt_A_buffer, AA[ass,:], bb[ass])
    # get_J_sd_wrt_b_fun(J_sd_wrt_b_buffer, AA[ass,:], bb[ass])

    get_sd_wrt_xt!(sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer)
    # get_J_sd_wrt_xt!(J_sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer, J_sd_wrt_A_buffer, J_A_wrt_xt_buffer, J_sd_wrt_b_buffer, J_b_wrt_xt_buffer)
    
    sd_lag .= [-λsd * sd_wrt_xt_buffer; zeros(6)]
end

# get_sd_lag!(sd_lag_buf, x0, λsd, i, ass, AA, bb)
get_sd_lag!(sd_lag_buf, x0, 1.1, i, ass, AA, bb)

# sd_lag = -λsd*∇sd, Jsdlag = [-λsd * Jsd -∇sd] ∈ R12 × R13
function get_Jsdlag!(Jsdlag, xt, λsd, i, ass, AA, bb)
    get_A_ego_wrt_xt_fun[i](A_ego_wrt_xt_buffer, xt[1:6])
    get_b_ego_wrt_xt_fun[i](b_ego_wrt_xt_buffer, xt[1:6])
    get_J_A_ego_wrt_xt_fun[i](J_A_ego_wrt_xt_buffer, xt[1:6])
    get_J_b_ego_wrt_xt_fun[i](J_b_ego_wrt_xt_buffer, xt[1:6])

    get_Ab_wrt_xt!(A_wrt_xt_buffer, b_wrt_xt_buffer, ass, A_ego_wrt_xt_buffer, b_ego_wrt_xt_buffer)
    get_J_Ab_wrt_xt!(J_A_wrt_xt_buffer, J_b_wrt_xt_buffer, ass, J_A_ego_wrt_xt_buffer, J_b_ego_wrt_xt_buffer)

    get_sd_wrt_A_fun(sd_wrt_A_buffer, AA[ass,:], bb[ass])
    get_sd_wrt_b_fun(sd_wrt_b_buffer, AA[ass,:], bb[ass])
    get_J_sd_wrt_A_fun(J_sd_wrt_A_buffer, AA[ass,:], bb[ass])
    get_J_sd_wrt_b_fun(J_sd_wrt_b_buffer, AA[ass,:], bb[ass])

    get_sd_wrt_xt!(sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer)
    get_J_sd_wrt_xt!(J_sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer, J_sd_wrt_A_buffer, J_A_wrt_xt_buffer, J_sd_wrt_b_buffer, J_b_wrt_xt_buffer)

    Jsdlag .= [
                [-λsd*J_sd_wrt_xt_buffer zeros(6, 6); zeros(6, 6)  zeros(6, 6)]     [-sd_wrt_xt_buffer; zeros(6)]
                ]
end
# get_Jsdlag!(Jsdlag, xt, λsd, i, ass, AA, bb)
get_Jsdlag!(Jsdlag_buf, x0, 1.1, i, ass, AA, bb)




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

