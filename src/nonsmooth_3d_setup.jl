
function get_Num0s(n_row, n_col)
    return [Num(0) for _ in 1:n_row, _ in 1:n_col]
end

function get_Num0s(n)
    return [Num(0) for _ in 1:n]
end

function get_aibi_wrt_xt_fun(R, l, xt)
    # Symbolics of one row of original ego constraints
    Ai = Symbolics.@variables(Ai[1:3])[1] |> Symbolics.scalarize
    bi = Symbolics.@variables(bi)[1]

    # parameterized Ai and bi w.r.t. xt
    a1, a2, a3 = Ai' * R'
    b1 = bi - [a1, a2, a3]' *l

    # expression of gradient
    grad_a1_wrt_xt = Symbolics.gradient(a1, xt; simplify=false)
    grad_a2_wrt_xt = Symbolics.gradient(a2, xt; simplify=false)
    grad_a3_wrt_xt = Symbolics.gradient(a3, xt; simplify=false)
    grad_b1_wrt_xt = Symbolics.gradient(b1, xt; simplify=false)

    Hessian_a1_wrt_xt = Symbolics.jacobian(grad_a1_wrt_xt, xt; simplify=false)
    Hessian_a2_wrt_xt = Symbolics.jacobian(grad_a2_wrt_xt, xt; simplify=false)
    Hessian_a3_wrt_xt = Symbolics.jacobian(grad_a3_wrt_xt, xt; simplify=false)
    Hessian_b1_wrt_xt = Symbolics.jacobian(grad_b1_wrt_xt, xt; simplify=false)

    # function of gradient
    get_grad_a1_wrt_xt = Symbolics.build_function(grad_a1_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    get_grad_a2_wrt_xt = Symbolics.build_function(grad_a2_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    get_grad_a3_wrt_xt = Symbolics.build_function(grad_a3_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    get_grad_b1_wrt_xt = Symbolics.build_function(grad_b1_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    
    # function of Hessian matrix
    get_Hessian_a1_wrt_xt = Symbolics.build_function(Hessian_a1_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    get_Hessian_a2_wrt_xt = Symbolics.build_function(Hessian_a2_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    get_Hessian_a3_wrt_xt = Symbolics.build_function(Hessian_a3_wrt_xt, xt, Ai, bi; expression=Val(false))[2]
    get_Hessian_b1_wrt_xt = Symbolics.build_function(Hessian_b1_wrt_xt, xt, Ai, bi; expression=Val(false))[2]

    return [get_grad_a1_wrt_xt, get_grad_a2_wrt_xt, get_grad_a3_wrt_xt, get_grad_b1_wrt_xt, get_Hessian_a1_wrt_xt, get_Hessian_a2_wrt_xt, get_Hessian_a3_wrt_xt, get_Hessian_b1_wrt_xt]
end

function get_Ab_ego_wrt_xt_fun(xt, A_ego, b_ego, aibi_wrt_xt_functions)
    get_grad_a1_wrt_xt, get_grad_a2_wrt_xt, get_grad_a3_wrt_xt, get_grad_b1_wrt_xt, get_Hessian_a1_wrt_xt, get_Hessian_a2_wrt_xt, get_Hessian_a3_wrt_xt, get_Hessian_b1_wrt_xt = aibi_wrt_xt_functions
    m1 = length(b_ego)

    grad_A_ego_wrt_xt = [get_Num0s(6) for _ in 1:m1, _ in 1:3] # the same shape as A_ego, where every element is a gradient
    grad_b_ego_wrt_xt = [get_Num0s(6) for _ in 1:m1] # the same shape as b_ego, where every element is a gradient 
    Hessian_A_ego_wrt_xt = [get_Num0s(6, 6) for _ in 1:m1, _ in 1:3] # the same shape as A_ego, where every element is a Hessian matrix
    Hessian_b_ego_wrt_xt = [get_Num0s(6, 6) for _ in 1:m1] # the same shape as b_ego, where every element is a Hessian matrix

    # replace symbolic Ai and bi with values of A_ego and b_ego, get expressions only including xt
    for k in 1:m1
        get_grad_a1_wrt_xt(grad_A_ego_wrt_xt[k, 1], xt, A_ego[k,:], b_ego[k])
        get_grad_a2_wrt_xt(grad_A_ego_wrt_xt[k, 2], xt, A_ego[k,:], b_ego[k])
        get_grad_a3_wrt_xt(grad_A_ego_wrt_xt[k, 3], xt, A_ego[k,:], b_ego[k])
        get_grad_b1_wrt_xt(grad_b_ego_wrt_xt[k], xt, A_ego[k,:], b_ego[k])

        get_Hessian_a1_wrt_xt(Hessian_A_ego_wrt_xt[k, 1], xt, A_ego[k,:], b_ego[k])
        get_Hessian_a2_wrt_xt(Hessian_A_ego_wrt_xt[k, 2], xt, A_ego[k,:], b_ego[k])
        get_Hessian_a3_wrt_xt(Hessian_A_ego_wrt_xt[k, 3], xt, A_ego[k,:], b_ego[k])
        get_Hessian_b1_wrt_xt(Hessian_b_ego_wrt_xt[k], xt, A_ego[k,:], b_ego[k])
    end

    get_grad_A_ego_wrt_xt = Symbolics.build_function(grad_A_ego_wrt_xt, xt; expression=Val(false))[2]
    get_grad_b_ego_wrt_xt = Symbolics.build_function(grad_b_ego_wrt_xt, xt; expression=Val(false))[2]
    get_Hessian_A_ego_wrt_xt = Symbolics.build_function(Hessian_A_ego_wrt_xt, xt; expression=Val(false))[2]
    get_Hessian_b_ego_wrt_xt = Symbolics.build_function(Hessian_b_ego_wrt_xt, xt; expression=Val(false))[2]
    return get_grad_A_ego_wrt_xt, get_grad_b_ego_wrt_xt, get_Hessian_A_ego_wrt_xt, get_Hessian_b_ego_wrt_xt
end

# update grad_A_wrt_xt and grad_b_wrt_xt according to assignment
function get_grad_Ab_wrt_xt!(grad_A_wrt_xt, grad_b_wrt_xt, ass, grad_A_ego_wrt_xt, grad_b_ego_wrt_xt, m1)
    for k in eachindex(grad_A_wrt_xt)
        grad_A_wrt_xt[k] .= .0
    end
    for k in eachindex(grad_b_wrt_xt)
        grad_b_wrt_xt[k] .= .0
    end

    for (k, ind) in enumerate(ass)
        if ind <= m1
            # must update element by element, otherwise they will point to the same memory, like the following line
            # A_wrt_xt[k, 1:3] .= A_ego_wrt_xt[ind, :]
            grad_A_wrt_xt[k, 1] .= grad_A_ego_wrt_xt[ind, 1]
            grad_A_wrt_xt[k, 2] .= grad_A_ego_wrt_xt[ind, 2]
            grad_A_wrt_xt[k, 3] .= grad_A_ego_wrt_xt[ind, 3]
            grad_b_wrt_xt[k] .= grad_b_ego_wrt_xt[ind]
        end
    end
end

# update Hessian_A_wrt_xt and Hessian_b_wrt_xt according to assignment
function get_Hessian_Ab_wrt_xt!(Hessian_A_wrt_xt, Hessian_b_wrt_xt, ass, Hessian_A_ego_wrt_xt, Hessian_b_ego_wrt_xt, m1)
    for k in eachindex(Hessian_A_wrt_xt)
        Hessian_A_wrt_xt[k] .= .0
    end
    for k in eachindex(Hessian_b_wrt_xt)
        Hessian_b_wrt_xt[k] .= .0
    end

    for (k, ind) in enumerate(ass)
        if ind <= m1
            # must update element by element, otherwise they will point to the same memory, like the following line
            # J_A_wrt_xt[k, 1:3] .= J_A_ego_wrt_xt[ind, :]
            Hessian_A_wrt_xt[k, 1] .= Hessian_A_ego_wrt_xt[ind, 1]
            Hessian_A_wrt_xt[k, 2] .= Hessian_A_ego_wrt_xt[ind, 2]
            Hessian_A_wrt_xt[k, 3] .= Hessian_A_ego_wrt_xt[ind, 3]
            Hessian_b_wrt_xt[k] .= Hessian_b_ego_wrt_xt[ind]
        end
    end
end

function get_Jacobian_Ab_wrt_xt(grad_A_wrt_xt, grad_b_wrt_xt)
    Jacobian_Ab_wrt_xt = zeros(20, 6)
    for k in eachindex(grad_A_wrt_xt)
        Jacobian_Ab_wrt_xt[k, :] .= grad_A_wrt_xt[k]
    end
    for k in eachindex(grad_b_wrt_xt)
        Jacobian_Ab_wrt_xt[k+16, :] .= grad_b_wrt_xt[k]
    end

    return Jacobian_Ab_wrt_xt    
end

function get_sd_wrt_Ab_fun()
    # Symbolics of a 4d linear system
    A = Symbolics.@variables(A[1:4, 1:4])[1] |> Symbolics.scalarize
    b = Symbolics.@variables(b[1:4])[1] |> Symbolics.scalarize

    # sol = -A \ b
    # for k in eachindex(sol)
    #     sol[k] = Symbolics.simplify(sol[k]) # need to be simplified, otherwise it returns NaN for some values of A and b
    # end

    # results of code above
    sol = [(A[1, 2]*A[2, 3]*A[3, 4]*b[4] - A[1, 2]*A[2, 3]*A[4, 4]*b[3] - A[1, 2]*A[2, 4]*A[3, 3]*b[4] + A[1, 2]*A[2, 4]*A[4, 3]*b[3] + A[1, 2]*A[3, 3]*A[4, 4]*b[2] - A[1, 2]*A[3, 4]*A[4, 3]*b[2] - A[1, 3]*A[2, 2]*A[3, 4]*b[4] + A[1, 3]*A[2, 2]*A[4, 4]*b[3] + A[1, 3]*A[2, 4]*A[3, 2]*b[4] - A[1, 3]*A[2, 4]*A[4, 2]*b[3] - A[1, 3]*A[3, 2]*A[4, 4]*b[2] + A[1, 3]*A[3, 4]*A[4, 2]*b[2] + A[1, 4]*A[2, 2]*A[3, 3]*b[4] - A[1, 4]*A[2, 2]*A[4, 3]*b[3] - A[1, 4]*A[2, 3]*A[3, 2]*b[4] + A[1, 4]*A[2, 3]*A[4, 2]*b[3] + A[1, 4]*A[3, 2]*A[4, 3]*b[2] - A[1, 4]*A[3, 3]*A[4, 2]*b[2] - A[2, 2]*A[3, 3]*A[4, 4]*b[1] + A[2, 2]*A[3, 4]*A[4, 3]*b[1] + A[2, 3]*A[3, 2]*A[4, 4]*b[1] - A[2, 3]*A[3, 4]*A[4, 2]*b[1] - A[2, 4]*A[3, 2]*A[4, 3]*b[1] + A[2, 4]*A[3, 3]*A[4, 2]*b[1]) / (A[1, 1]*A[2, 2]*A[3, 3]*A[4, 4] - A[1, 1]*A[2, 2]*A[3, 4]*A[4, 3] - A[1, 1]*A[2, 3]*A[3, 2]*A[4, 4] + A[1, 1]*A[2, 3]*A[3, 4]*A[4, 2] + A[1, 1]*A[2, 4]*A[3, 2]*A[4, 3] - A[1, 1]*A[2, 4]*A[3, 3]*A[4, 2] - A[1, 2]*A[2, 1]*A[3, 3]*A[4, 4] + A[1, 2]*A[2, 1]*A[3, 4]*A[4, 3] + A[1, 2]*A[2, 3]*A[3, 1]*A[4, 4] - A[1, 2]*A[2, 3]*A[3, 4]*A[4, 1] - A[1, 2]*A[2, 4]*A[3, 1]*A[4, 3] + A[1, 2]*A[2, 4]*A[3, 3]*A[4, 1] + A[1, 3]*A[2, 1]*A[3, 2]*A[4, 4] - A[1, 3]*A[2, 1]*A[3, 4]*A[4, 2] - A[1, 3]*A[2, 2]*A[3, 1]*A[4, 4] + A[1, 3]*A[2, 2]*A[3, 4]*A[4, 1] + A[1, 3]*A[2, 4]*A[3, 1]*A[4, 2] - A[1, 3]*A[2, 4]*A[3, 2]*A[4, 1] - A[1, 4]*A[2, 1]*A[3, 2]*A[4, 3] + A[1, 4]*A[2, 1]*A[3, 3]*A[4, 2] + A[1, 4]*A[2, 2]*A[3, 1]*A[4, 3] - A[1, 4]*A[2, 2]*A[3, 3]*A[4, 1] - A[1, 4]*A[2, 3]*A[3, 1]*A[4, 2] + A[1, 4]*A[2, 3]*A[3, 2]*A[4, 1])
            (-A[1, 1]*A[2, 3]*A[3, 4]*b[4] + A[1, 1]*A[2, 3]*A[4, 4]*b[3] + A[1, 1]*A[2, 4]*A[3, 3]*b[4] - A[1, 1]*A[2, 4]*A[4, 3]*b[3] - A[1, 1]*A[3, 3]*A[4, 4]*b[2] + A[1, 1]*A[3, 4]*A[4, 3]*b[2] + A[1, 3]*A[2, 1]*A[3, 4]*b[4] - A[1, 3]*A[2, 1]*A[4, 4]*b[3] - A[1, 3]*A[2, 4]*A[3, 1]*b[4] + A[1, 3]*A[2, 4]*A[4, 1]*b[3] + A[1, 3]*A[3, 1]*A[4, 4]*b[2] - A[1, 3]*A[3, 4]*A[4, 1]*b[2] - A[1, 4]*A[2, 1]*A[3, 3]*b[4] + A[1, 4]*A[2, 1]*A[4, 3]*b[3] + A[1, 4]*A[2, 3]*A[3, 1]*b[4] - A[1, 4]*A[2, 3]*A[4, 1]*b[3] - A[1, 4]*A[3, 1]*A[4, 3]*b[2] + A[1, 4]*A[3, 3]*A[4, 1]*b[2] + A[2, 1]*A[3, 3]*A[4, 4]*b[1] - A[2, 1]*A[3, 4]*A[4, 3]*b[1] - A[2, 3]*A[3, 1]*A[4, 4]*b[1] + A[2, 3]*A[3, 4]*A[4, 1]*b[1] + A[2, 4]*A[3, 1]*A[4, 3]*b[1] - A[2, 4]*A[3, 3]*A[4, 1]*b[1]) / (A[1, 1]*A[2, 2]*A[3, 3]*A[4, 4] - A[1, 1]*A[2, 2]*A[3, 4]*A[4, 3] - A[1, 1]*A[2, 3]*A[3, 2]*A[4, 4] + A[1, 1]*A[2, 3]*A[3, 4]*A[4, 2] + A[1, 1]*A[2, 4]*A[3, 2]*A[4, 3] - A[1, 1]*A[2, 4]*A[3, 3]*A[4, 2] - A[1, 2]*A[2, 1]*A[3, 3]*A[4, 4] + A[1, 2]*A[2, 1]*A[3, 4]*A[4, 3] + A[1, 2]*A[2, 3]*A[3, 1]*A[4, 4] - A[1, 2]*A[2, 3]*A[3, 4]*A[4, 1] - A[1, 2]*A[2, 4]*A[3, 1]*A[4, 3] + A[1, 2]*A[2, 4]*A[3, 3]*A[4, 1] + A[1, 3]*A[2, 1]*A[3, 2]*A[4, 4] - A[1, 3]*A[2, 1]*A[3, 4]*A[4, 2] - A[1, 3]*A[2, 2]*A[3, 1]*A[4, 4] + A[1, 3]*A[2, 2]*A[3, 4]*A[4, 1] + A[1, 3]*A[2, 4]*A[3, 1]*A[4, 2] - A[1, 3]*A[2, 4]*A[3, 2]*A[4, 1] - A[1, 4]*A[2, 1]*A[3, 2]*A[4, 3] + A[1, 4]*A[2, 1]*A[3, 3]*A[4, 2] + A[1, 4]*A[2, 2]*A[3, 1]*A[4, 3] - A[1, 4]*A[2, 2]*A[3, 3]*A[4, 1] - A[1, 4]*A[2, 3]*A[3, 1]*A[4, 2] + A[1, 4]*A[2, 3]*A[3, 2]*A[4, 1])
            (A[1, 1]*A[2, 2]*A[3, 4]*b[4] - A[1, 1]*A[2, 2]*A[4, 4]*b[3] - A[1, 1]*A[2, 4]*A[3, 2]*b[4] + A[1, 1]*A[2, 4]*A[4, 2]*b[3] + A[1, 1]*A[3, 2]*A[4, 4]*b[2] - A[1, 1]*A[3, 4]*A[4, 2]*b[2] - A[1, 2]*A[2, 1]*A[3, 4]*b[4] + A[1, 2]*A[2, 1]*A[4, 4]*b[3] + A[1, 2]*A[2, 4]*A[3, 1]*b[4] - A[1, 2]*A[2, 4]*A[4, 1]*b[3] - A[1, 2]*A[3, 1]*A[4, 4]*b[2] + A[1, 2]*A[3, 4]*A[4, 1]*b[2] + A[1, 4]*A[2, 1]*A[3, 2]*b[4] - A[1, 4]*A[2, 1]*A[4, 2]*b[3] - A[1, 4]*A[2, 2]*A[3, 1]*b[4] + A[1, 4]*A[2, 2]*A[4, 1]*b[3] + A[1, 4]*A[3, 1]*A[4, 2]*b[2] - A[1, 4]*A[3, 2]*A[4, 1]*b[2] - A[2, 1]*A[3, 2]*A[4, 4]*b[1] + A[2, 1]*A[3, 4]*A[4, 2]*b[1] + A[2, 2]*A[3, 1]*A[4, 4]*b[1] - A[2, 2]*A[3, 4]*A[4, 1]*b[1] - A[2, 4]*A[3, 1]*A[4, 2]*b[1] + A[2, 4]*A[3, 2]*A[4, 1]*b[1]) / (A[1, 1]*A[2, 2]*A[3, 3]*A[4, 4] - A[1, 1]*A[2, 2]*A[3, 4]*A[4, 3] - A[1, 1]*A[2, 3]*A[3, 2]*A[4, 4] + A[1, 1]*A[2, 3]*A[3, 4]*A[4, 2] + A[1, 1]*A[2, 4]*A[3, 2]*A[4, 3] - A[1, 1]*A[2, 4]*A[3, 3]*A[4, 2] - A[1, 2]*A[2, 1]*A[3, 3]*A[4, 4] + A[1, 2]*A[2, 1]*A[3, 4]*A[4, 3] + A[1, 2]*A[2, 3]*A[3, 1]*A[4, 4] - A[1, 2]*A[2, 3]*A[3, 4]*A[4, 1] - A[1, 2]*A[2, 4]*A[3, 1]*A[4, 3] + A[1, 2]*A[2, 4]*A[3, 3]*A[4, 1] + A[1, 3]*A[2, 1]*A[3, 2]*A[4, 4] - A[1, 3]*A[2, 1]*A[3, 4]*A[4, 2] - A[1, 3]*A[2, 2]*A[3, 1]*A[4, 4] + A[1, 3]*A[2, 2]*A[3, 4]*A[4, 1] + A[1, 3]*A[2, 4]*A[3, 1]*A[4, 2] - A[1, 3]*A[2, 4]*A[3, 2]*A[4, 1] - A[1, 4]*A[2, 1]*A[3, 2]*A[4, 3] + A[1, 4]*A[2, 1]*A[3, 3]*A[4, 2] + A[1, 4]*A[2, 2]*A[3, 1]*A[4, 3] - A[1, 4]*A[2, 2]*A[3, 3]*A[4, 1] - A[1, 4]*A[2, 3]*A[3, 1]*A[4, 2] + A[1, 4]*A[2, 3]*A[3, 2]*A[4, 1])
            (-A[1, 1]*A[2, 2]*A[3, 3]*b[4] + A[1, 1]*A[2, 2]*A[4, 3]*b[3] + A[1, 1]*A[2, 3]*A[3, 2]*b[4] - A[1, 1]*A[2, 3]*A[4, 2]*b[3] - A[1, 1]*A[3, 2]*A[4, 3]*b[2] + A[1, 1]*A[3, 3]*A[4, 2]*b[2] + A[1, 2]*A[2, 1]*A[3, 3]*b[4] - A[1, 2]*A[2, 1]*A[4, 3]*b[3] - A[1, 2]*A[2, 3]*A[3, 1]*b[4] + A[1, 2]*A[2, 3]*A[4, 1]*b[3] + A[1, 2]*A[3, 1]*A[4, 3]*b[2] - A[1, 2]*A[3, 3]*A[4, 1]*b[2] - A[1, 3]*A[2, 1]*A[3, 2]*b[4] + A[1, 3]*A[2, 1]*A[4, 2]*b[3] + A[1, 3]*A[2, 2]*A[3, 1]*b[4] - A[1, 3]*A[2, 2]*A[4, 1]*b[3] - A[1, 3]*A[3, 1]*A[4, 2]*b[2] + A[1, 3]*A[3, 2]*A[4, 1]*b[2] + A[2, 1]*A[3, 2]*A[4, 3]*b[1] - A[2, 1]*A[3, 3]*A[4, 2]*b[1] - A[2, 2]*A[3, 1]*A[4, 3]*b[1] + A[2, 2]*A[3, 3]*A[4, 1]*b[1] + A[2, 3]*A[3, 1]*A[4, 2]*b[1] - A[2, 3]*A[3, 2]*A[4, 1]*b[1]) / (A[1, 1]*A[2, 2]*A[3, 3]*A[4, 4] - A[1, 1]*A[2, 2]*A[3, 4]*A[4, 3] - A[1, 1]*A[2, 3]*A[3, 2]*A[4, 4] + A[1, 1]*A[2, 3]*A[3, 4]*A[4, 2] + A[1, 1]*A[2, 4]*A[3, 2]*A[4, 3] - A[1, 1]*A[2, 4]*A[3, 3]*A[4, 2] - A[1, 2]*A[2, 1]*A[3, 3]*A[4, 4] + A[1, 2]*A[2, 1]*A[3, 4]*A[4, 3] + A[1, 2]*A[2, 3]*A[3, 1]*A[4, 4] - A[1, 2]*A[2, 3]*A[3, 4]*A[4, 1] - A[1, 2]*A[2, 4]*A[3, 1]*A[4, 3] + A[1, 2]*A[2, 4]*A[3, 3]*A[4, 1] + A[1, 3]*A[2, 1]*A[3, 2]*A[4, 4] - A[1, 3]*A[2, 1]*A[3, 4]*A[4, 2] - A[1, 3]*A[2, 2]*A[3, 1]*A[4, 4] + A[1, 3]*A[2, 2]*A[3, 4]*A[4, 1] + A[1, 3]*A[2, 4]*A[3, 1]*A[4, 2] - A[1, 3]*A[2, 4]*A[3, 2]*A[4, 1] - A[1, 4]*A[2, 1]*A[3, 2]*A[4, 3] + A[1, 4]*A[2, 1]*A[3, 3]*A[4, 2] + A[1, 4]*A[2, 2]*A[3, 1]*A[4, 3] - A[1, 4]*A[2, 2]*A[3, 3]*A[4, 1] - A[1, 4]*A[2, 3]*A[3, 1]*A[4, 2] + A[1, 4]*A[2, 3]*A[3, 2]*A[4, 1])]
    sd = sol[4]
    # symbolic solution to a 4d linear system
    get_sol_A_b = Symbolics.build_function(sol, A, b; expression=Val(false))[2]

    Ab = vcat(vec(A),b) # vectorize all parameters
    grad_sd_wrt_Ab = Symbolics.gradient(sd, Ab)
    Hessian_sd_wrt_Ab = Symbolics.jacobian(grad_sd_wrt_Ab, Ab)
    get_grad_sd_wrt_Ab = Symbolics.build_function(grad_sd_wrt_Ab, Ab; expression=Val(false))[2]
    get_Hessian_sd_wrt_Ab = Symbolics.build_function(Hessian_sd_wrt_Ab, Ab; expression=Val(false))[2]

    return get_grad_sd_wrt_Ab, get_Hessian_sd_wrt_Ab, get_sol_A_b
end

function get_sol_xt!(sol_xt, AA, bb, get_sol_A_b, assignments)
    sol = zeros(4)
    for (k, ass) in enumerate(assignments)
        get_sol_A_b(sol, AA[ass, :], bb[ass])
        sol_xt[:, k] .= sol
    end
end

function get_grad_sd_wrt_xt!(grad_sd_wrt_xt, grad_sd_wrt_Ab, grad_A_wrt_xt, grad_b_wrt_xt)
    grad_sd_wrt_xt .= grad_sd_wrt_Ab' * vec([grad_A_wrt_xt grad_b_wrt_xt])
end

function get_Hessian_sd_wrt_xt!(Hessian_sd_wrt_xt, grad_sd_wrt_Ab, grad_A_wrt_xt, grad_b_wrt_xt, Hessian_sd_wrt_Ab, Hessian_A_wrt_xt, Hessian_b_wrt_xt)
    Hessian_sd_wrt_xt .= 0.0
    Jacobian_Ab_wrt_xt = get_Jacobian_Ab_wrt_xt(grad_A_wrt_xt, grad_b_wrt_xt)
    HJ = Hessian_sd_wrt_Ab * Jacobian_Ab_wrt_xt
    for k in eachindex(grad_sd_wrt_Ab)
        Hessian_sd_wrt_xt .+= Jacobian_Ab_wrt_xt[k, :] * HJ[k, :]' + grad_sd_wrt_Ab[k] * vec([Hessian_A_wrt_xt Hessian_b_wrt_xt])[k]
    end
end

function gen_LP_data_3d(A_ego::AbstractArray{T}, b_ego, centr_ego, A_obs, b_obs, centr_obs) where {T}
    A = [A_ego A_ego*centr_ego+b_ego
        A_obs A_obs*centr_obs+b_obs]
    b = [b_ego; b_obs]
    q = [0, 0, 0, 1.0]
    (A, b, q)
end

# filter indices which are impossible to be active at the same time for one poly
function get_3_possible_constraint_ids(A, b; tol=1e-4)
    AA = Matrix(A)
    dim = size(AA)[2]
    bb = b
    ind = collect(1:length(bb))
    inds = powerset(ind) |> collect
    itr = [i for i in inds if length(i)==dim]
    feasible_inds=[]
    for i in itr
        try
            xx = - AA[i,:] \ bb[i]
            if all(AA*xx +bb .> -tol)
                push!(feasible_inds, i)
            end
        catch err
            if err isa LinearAlgebra.SingularException
                continue
            else
                # @warn(err)
            end
        end
    end
    feasible_inds
end

# get possible pair of indices, which means edges
function get_2_possible_constraint_ids(A, b, V; tol=1e-4)
    matrix_v_face = map(V) do v
        A * v + b .< tol
    end
    feasible_inds = []
    for i in eachindex(matrix_v_face)
        for j in i+1:length(matrix_v_face)
            common_faces = matrix_v_face[i] + matrix_v_face[j] .== 2
            if sum(common_faces) == 2
                push!(feasible_inds, findall(common_faces))
            end
        end
    end
    feasible_inds
end

# get_possible_assignments_3d(Pe.A, Pe.b, Pe.V, Po.A, Po.b, Po.V)
# enumerate possible assignments (2 indices from one poly, and 2 indices from the other; or 3 + 1; or 1 + 3)
function get_possible_assignments_3d(Ae, be, Ve, Ao, bo, Vo)
    m1 = length(be)
    m2 = length(bo)
    inds_e = get_3_possible_constraint_ids(Ae, be)
    inds_o = get_3_possible_constraint_ids(Ao, bo)
    for i in eachindex(inds_o)
        inds_o[i] += [m1, m1, m1]
    end
    Itr=[]
    for i in 1:m1
        for ind in inds_o
            push!(Itr, sort(vcat(ind, i)))
        end
    end
    for i in m1+1:m1+m2
        for ind in inds_e
            push!(Itr, sort(vcat(ind, i)))
        end
    end

    inds_e = get_2_possible_constraint_ids(Ae, be, Ve)
    inds_o = get_2_possible_constraint_ids(Ao, bo, Vo)
    for i in eachindex(inds_o)
        inds_o[i] += [m1, m1]
    end
    for e in inds_e
        for o in inds_o
            push!(Itr, [e; o])
        end
    end

    Itr
end

function g_col_single_3d(xt, A_ego, b_ego, V_ego, centr_ego, A_obs, b_obs, V_obs, centr_obs)
    sds = Dict()
    intercepts = Dict()
    R = R_from_mrp(xt[4:6])
    l = xt[1:3]

    # some terms are not eliminated when calculating symbolics
    # Aex, bex = shift_to_3D(A_ego, b_ego, xt)
    # centroidex = l + R * centr_ego
    # AA, bb, qq = gen_LP_data_3d(Aex, bex, centroidex, A_obs, b_obs, centr_obs)

    # better way to get AA and bb than above
    AA = [A_ego*R' A_ego*centr_ego+b_ego
          A_obs    A_obs*centr_obs+b_obs]
    bb = [b_ego-A_ego*R'*l
          b_obs]
    Itr = get_possible_assignments_3d(A_ego, b_ego, V_ego, A_obs, b_obs, V_obs)
    for active_inds in Itr
        # if !is_ass_feasible(active_inds, m1, m2)
        #     continue
        # end

        try
            AA_active = collect(AA[active_inds, :])
            bb_active = collect(bb[active_inds])
            # TODO what if not unique primal? Need to resolve
            if length(active_inds) == 4
                zz = -AA_active \ bb_active
            else
                # if linear system is underdetermined, use minimum norm solution (calculated by right inverse)
                # Note: every solution has the same sd, i.e., zz[3], but different zz[1:2]
                zz = -AA_active' * ((AA_active * AA_active') \ bb_active)
                @warn "wrong"
                @infiltrate
            end
            sd = zz[4]
            intercepts[active_inds] = zz[1:3]
            sds[active_inds] = sd
        catch err
            if err isa LinearAlgebra.SingularException
                continue
            else
                @warn(err)
            end
        end
    end
    sds, intercepts, AA, bb
end

function setup_nonsmooth_3d(
    ego_polys,
    obs_polys;
    T=2,
    dt=0.2,
    R_cost=1e-3 * I(6), # penality for control variables
    Q_cost=1e-3 * I(3), # penality for distance
    p1_max=500.0,
    p2_max=500.0,
    p3_max=500.0,
    u1_max=1.0,
    u2_max=1.0,
    u3_max=1.0,
    u4_max=π / 4,
    u5_max=π / 4,
    u6_max=π / 4,
    n_sd_slots=2
)

    # problem dimensions
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
    n_sd_cons = T * n_ego * n_obs * n_sd_slots

    # assume number of sides are the same for all ego polys, obs polys
    # @assert all([length(ego_polys[i].b) for i in 1:n_ego] .== n_side_ego)
    # @assert all([length(obs_polys[i].b) for i in 1:n_obs] .== n_side_obs)

    # θ indexing
    z_s2i, dyn_cons_s2i, env_cons_s2i, sd_cons_s2i = get_sub2idxs((n_xu, T), (n_dyn_cons), (n_env_cons), (n_sd_slots, n_obs, n_ego, T))

    z = Symbolics.@variables(z[1:n_z])[1] |> Symbolics.scalarize
    x0 = Symbolics.@variables(x0[1:n_x])[1] |> Symbolics.scalarize
    xt = Symbolics.@variables(xt[1:n_x])[1] |> Symbolics.scalarize

    cost = f_3d(z, T, R_cost, Q_cost)
    dyn_cons = g_dyn_3d(z, x0, T, dt)
    env_cons = g_env_3d(z, T, p1_max, p2_max, p3_max, u1_max, u2_max, u3_max, u4_max, u5_max, u6_max)

    # check indexing consistency
    @assert length(z_s2i) == n_z
    @assert length(dyn_cons_s2i) == length(dyn_cons)
    @assert length(env_cons_s2i) == length(env_cons)
    @assert length(sd_cons_s2i) == n_sd_cons

    # sds, intercepts, AAs, bbs, etc to fill F and J
    get_AA = Dict()
    get_bb = Dict()
    Itr_dict = Dict()

    # functions vary according to different egos
    get_grad_A_ego_wrt_xt = Dict()
    get_grad_b_ego_wrt_xt = Dict()
    get_Hessian_A_ego_wrt_xt = Dict()
    get_Hessian_b_ego_wrt_xt = Dict()
    # we solve sds symbolically for given ego and obs at problem creation, 
    # TODO could be done in a more flexible by abstracting problem parameters (obs_polys) and filling them in fill_F, fill_J instead as we did before
    λsd = Symbolics.@variables(λsd)[1]

    # rotation matrix and translation 
    R = R_from_mrp(xt[4:6])
    l = xt[1:3]

    # helper functions to calculate Jacobian
    aibi_wrt_xt_functions = get_aibi_wrt_xt_fun(R, l, xt[1:6])
    get_grad_sd_wrt_Ab, get_Hessian_sd_wrt_Ab, get_sol_A_b = get_sd_wrt_Ab_fun()

    #@info "Generating symbolic sds, intercepts, AAs, bbs"
    #p = Progress(n_sds * n_ego * n_obs, dt=1.0)
    for (i, Pe) in enumerate(ego_polys)
        A_ego = collect(Pe.A)
        b_ego = Pe.b
        V_ego = Pe.V
        centr_ego = Pe.c
        get_grad_A_ego_wrt_xt[i], get_grad_b_ego_wrt_xt[i], get_Hessian_A_ego_wrt_xt[i], get_Hessian_b_ego_wrt_xt[i] = get_Ab_ego_wrt_xt_fun(xt[1:6], A_ego, b_ego, aibi_wrt_xt_functions)

        for (j, Po) in enumerate(obs_polys)
            A_obs = collect(Po.A)
            b_obs = Po.b
            V_obs = Po.V
            centr_obs = Po.c
            
            AA = [A_ego*R' A_ego*centr_ego+b_ego
                A_obs    A_obs*centr_obs+b_obs]
            bb = [b_ego-A_ego*R'*l
                b_obs]

            Itr_dict[i, j] = get_possible_assignments_3d(A_ego, b_ego, V_ego, A_obs, b_obs, V_obs)
            get_AA[i, j] = Symbolics.build_function(AA', xt; expression=Val(false))[2]
            get_bb[i, j] = Symbolics.build_function(bb, xt; expression=Val(false))[2]

        end
    end
    # different ego-obs pair may have different numbers of possible assignments
    n_sds = []
    for (i, dic) in enumerate(Itr_dict)
        push!(n_sds, length(dic[2]))
    end
    n_sds = maximum(n_sds)

    nom_cons = [dyn_cons; env_cons]
    λ_nom = Symbolics.@variables(λ_nom[1:length(nom_cons)])[1] |> Symbolics.scalarize
    λ_sd = Symbolics.@variables(λ_col[1:n_sd_cons])[1] |> Symbolics.scalarize

    # F = [F_nom; F_sd (to be filled later)]
    # J = dJ/dθ
    # we also need sds/dx
# @infiltrate
    # F and J nominal (before filling)
    θ = [z; λ_nom; λ_sd]
    lag = cost - nom_cons' * λ_nom #- λ_sd' * nom_sd (to be filled later)
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; nom_cons; zeros(Num, n_sd_cons)]
    J_nom = Symbolics.sparsejacobian(F_nom, θ)
# @infiltrate
    n = length(F_nom)
    l = zeros(n)
    u = zeros(n)

    # gradient = 0
    l[z_s2i[:]] .= -Inf
    u[z_s2i[:]] .= Inf

    # dynamics constraints
    l[dyn_cons_s2i[:]] .= -Inf
    u[dyn_cons_s2i[:]] .= Inf

    # environmental constraints
    l[env_cons_s2i[:]] .= 0.0
    u[env_cons_s2i[:]] .= Inf

    # sd constraints
    l[sd_cons_s2i[:]] .= 0.0
    u[sd_cons_s2i[:]] .= Inf

    # get nominal F and J (later to be filled by fill_F! and fill_J!)
    get_Fnom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false))[2]
# @infiltrate
    (Jnom_rows, Jnom_cols, Jnom_vals) = findnz(J_nom)
    get_Jnom_vals! = Symbolics.build_function(Jnom_vals, z, x0, λ_nom; expression=Val(false))[2]

    # sort sds and intercepts for given xt
    sol_xt_buffer_full = zeros(4, n_sds)
    AA_buffer_full = zeros(4, n_side_ego + n_side_obs)
    bb_buffer_full = zeros(n_side_ego + n_side_obs)

    function get_sorted_sds_3d(i, m1, j, m2, xt; tol=1e-4, local_factor=1.5)

        # if the size of buffer is larger than what the function returns, it is filled by column, so here AA is transposed
        get_AA[i, j](AA_buffer_full, xt)
        get_bb[i, j](bb_buffer_full, xt)
        AA_buffer = AA_buffer_full'[1:m1+m2, :]
        bb_buffer = bb_buffer_full[1:m1+m2]

        assignments = Itr_dict[i, j]
        n_ass = length(assignments)
        get_sol_xt!(sol_xt_buffer_full, AA_buffer, bb_buffer, get_sol_A_b, assignments)
        sol_xt_buffer = @view sol_xt_buffer_full'[1:n_ass, :]
        
        sds_buffer = @view sol_xt_buffer[:, end]
        intercept_buffer = @view sol_xt_buffer[:, 1:end-1]
# @infiltrate
        # tol = 1e-4
        # sd must be greater than -1 to be valid
        valid_mask = sds_buffer .>= -1.0 - tol

        # sd and intercept must correspond to a feasible vertex to be valid
        zz_check = hcat(intercept_buffer[valid_mask, :], sds_buffer[valid_mask])

        valid_mask[valid_mask] = map(eachrow(zz_check)) do row
            all(AA_buffer * row + bb_buffer .>= -tol)
        end

        sorted_sds_inds = sds_buffer[valid_mask] |> sortperm
        sorted_sds = sds_buffer[valid_mask][sorted_sds_inds]
        sorted_ass = assignments[valid_mask][sorted_sds_inds]

        # need to consider smarter way to filter out distant sds
        # local_factor = 1.5 # regard sds which are less than sd*local_factor as potential true sds
        local_sd_mask = (sorted_sds .+ 1) .<= (sorted_sds[1] + 1) * local_factor
        sorted_sds = sorted_sds[local_sd_mask]
        sorted_ass = sorted_ass[local_sd_mask]
        (sorted_sds, sorted_ass, AA_buffer, bb_buffer)
    end

    # compute the most consistent mapping between sd slot memory and current assignments
    sd_slot_mem = [Vector{Int64}() for _ in 1:n_sd_slots]

    function compute_ass_ind_map(sorted_ass, n_ass)
        ind_dict = Dict()
        ind_map = collect(1:n_ass)

        # identify which assignments exist in sd_slot_mem
        #for i in 1:n_ass
        #    for (k, mem) in enumerate(sd_slot_mem)
        #        if sorted_ass[i] == mem
        #            ind_dict[i] = k
        #        end
        #    end
        #end

        ## place remaining keys arbitrarily
        #for i in 1:n_ass
        #    if !haskey(ind_dict, i)
        #        k = 1
        #        while k ∈ values(ind_dict)
        #            k += 1
        #        end

        #        ind_dict[i] = k
        #    end
        #end

        #for (k, i) in ind_dict
        #    ind_map[k] = i
        #end

        #@infiltrate any(ind_map .!= collect(1:n_ass))
        ind_map = collect(1:n_ass) # smart indexing disabled
        ind_map
    end

    # fill_F!
    λ_nom_s2i = [dyn_cons_s2i...; env_cons_s2i...]

        

    # do not use fill(), because every element points to the same vector zeros(6)
    grad_A_ego_wrt_xt_buffer = [zeros(6) for _ in 1:n_side_ego, _ in 1:3]
    grad_b_ego_wrt_xt_buffer = [zeros(6) for _ in 1:n_side_ego]
    Hessian_A_ego_wrt_xt_buffer = [zeros(6, 6) for _ in 1:n_side_ego, _ in 1:3]
    Hessian_b_ego_wrt_xt_buffer = [zeros(6, 6) for _ in 1:n_side_ego]

    grad_A_wrt_xt_buffer = [zeros(6) for _ in 1:4, _ in 1:4]
    grad_b_wrt_xt_buffer = [zeros(6) for _ in 1:4]
    Hessian_A_wrt_xt_buffer = [zeros(6, 6) for _ in 1:4, _ in 1:4]
    Hessian_b_wrt_xt_buffer = [zeros(6, 6) for _ in 1:4]

    grad_sd_wrt_Ab_buffer = zeros(20)
    Hessian_sd_wrt_Ab_buffer = zeros(20, 20)

    # grad_sd_wrt_A_buffer = zeros(4,4)
    # grad_sd_wrt_b_buffer = zeros(4)
    # Hessian_sd_wrt_A_buffer = zeros(4,4)
    # Hessian_sd_wrt_b_buffer = zeros(4)

    grad_sd_wrt_xt_buffer = zeros(6)
    Hessian_sd_wrt_xt_buffer = zeros(6, 6)

    sd_lag_buf = zeros(n_x)
    Jsdlag_buf = zeros(n_x, n_x+1)
    Jsd_buf = zeros(n_x)


    # # Jsd = ∇sd, 12d vector
    # function get_Jsd!(J_sd, xt, i, ass, AA, bb, m1)
    #     get_grad_A_ego_wrt_xt[i](grad_A_ego_wrt_xt_buffer, xt[1:6])
    #     get_grad_b_ego_wrt_xt[i](grad_b_ego_wrt_xt_buffer, xt[1:6])
    #     # get_J_A_ego_wrt_xt_fun[i](J_A_ego_wrt_xt_buffer, xt[1:6])
    #     # get_J_b_ego_wrt_xt_fun[i](J_b_ego_wrt_xt_buffer, xt[1:6])

    #     get_grad_Ab_wrt_xt!(grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer, ass, grad_A_ego_wrt_xt_buffer, grad_b_ego_wrt_xt_buffer, m1)
    #     # get_J_Ab_wrt_xt!(J_A_wrt_xt_buffer, J_b_wrt_xt_buffer, ass, J_A_ego_wrt_xt_buffer, J_b_ego_wrt_xt_buffer, m1)

    #     get_sd_wrt_A_fun(grad_sd_wrt_A_buffer, AA[ass,:], bb[ass])
    #     get_sd_wrt_b_fun(grad_sd_wrt_b_buffer, AA[ass,:], bb[ass])
    #     # get_J_sd_wrt_A_fun(J_sd_wrt_A_buffer, AA[ass,:], bb[ass])
    #     # get_J_sd_wrt_b_fun(J_sd_wrt_b_buffer, AA[ass,:], bb[ass])

    #     get_grad_sd_wrt_xt!(grad_sd_wrt_xt_buffer, grad_sd_wrt_A_buffer, grad_A_wrt_xt_buffer, grad_sd_wrt_b_buffer, grad_b_wrt_xt_buffer)
    #     # get_J_sd_wrt_xt!(J_sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer, J_sd_wrt_A_buffer, J_A_wrt_xt_buffer, J_sd_wrt_b_buffer, J_b_wrt_xt_buffer)
        
    #     J_sd .= [grad_sd_wrt_xt_buffer; zeros(6)]
    # end


    # # sd_lag = -λsd*∇sd, 12d vector
    # function get_sd_lag!(sd_lag, xt, λsd, i, ass, AA, bb, m1)
    #     get_grad_A_ego_wrt_xt[i](grad_A_ego_wrt_xt_buffer, xt[1:6])
    #     get_grad_b_ego_wrt_xt[i](grad_b_ego_wrt_xt_buffer, xt[1:6])
    #     # get_J_A_ego_wrt_xt_fun[i](J_A_ego_wrt_xt_buffer, xt[1:6])
    #     # get_J_b_ego_wrt_xt_fun[i](J_b_ego_wrt_xt_buffer, xt[1:6])

    #     get_grad_Ab_wrt_xt!(grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer, ass, grad_A_ego_wrt_xt_buffer, grad_b_ego_wrt_xt_buffer, m1)
    #     # get_J_Ab_wrt_xt!(J_A_wrt_xt_buffer, J_b_wrt_xt_buffer, ass, J_A_ego_wrt_xt_buffer, J_b_ego_wrt_xt_buffer, m1)

    #     get_sd_wrt_A_fun(grad_sd_wrt_A_buffer, AA[ass,:], bb[ass])
    #     get_sd_wrt_b_fun(grad_sd_wrt_b_buffer, AA[ass,:], bb[ass])
    #     # get_J_sd_wrt_A_fun(J_sd_wrt_A_buffer, AA[ass,:], bb[ass])
    #     # get_J_sd_wrt_b_fun(J_sd_wrt_b_buffer, AA[ass,:], bb[ass])

    #     get_grad_sd_wrt_xt!(grad_sd_wrt_xt_buffer, grad_sd_wrt_A_buffer, grad_A_wrt_xt_buffer, grad_sd_wrt_b_buffer, grad_b_wrt_xt_buffer)
    #     # get_J_sd_wrt_xt!(J_sd_wrt_xt_buffer, sd_wrt_A_buffer, A_wrt_xt_buffer, sd_wrt_b_buffer, b_wrt_xt_buffer, J_sd_wrt_A_buffer, J_A_wrt_xt_buffer, J_sd_wrt_b_buffer, J_b_wrt_xt_buffer)
        
    #     sd_lag .= [-λsd * grad_sd_wrt_xt_buffer; zeros(6)]
    # end


    # # sd_lag = -λsd*∇sd, Jsdlag = [-λsd * Jsd -∇sd] ∈ R12 × R13
    # function get_Jsdlag!(Jsdlag, xt, λsd, i, ass, AA, bb, m1)
    #     get_grad_A_ego_wrt_xt[i](grad_A_ego_wrt_xt_buffer, xt[1:6])
    #     get_grad_b_ego_wrt_xt[i](grad_b_ego_wrt_xt_buffer, xt[1:6])
    #     get_Hessian_A_ego_wrt_xt[i](Hessian_A_ego_wrt_xt_buffer, xt[1:6])
    #     get_Hessian_b_ego_wrt_xt[i](Hessian_b_ego_wrt_xt_buffer, xt[1:6])

    #     get_grad_Ab_wrt_xt!(grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer, ass, grad_A_ego_wrt_xt_buffer, grad_b_ego_wrt_xt_buffer, m1)
    #     get_Hessian_Ab_wrt_xt!(Hessian_A_wrt_xt_buffer, Hessian_b_wrt_xt_buffer, ass, Hessian_A_ego_wrt_xt_buffer, Hessian_b_ego_wrt_xt_buffer, m1)

    #     get_sd_wrt_A_fun(grad_sd_wrt_A_buffer, AA[ass,:], bb[ass])
    #     get_sd_wrt_b_fun(grad_sd_wrt_b_buffer, AA[ass,:], bb[ass])
    #     get_J_sd_wrt_A_fun(Hessian_sd_wrt_A_buffer, AA[ass,:], bb[ass])
    #     get_J_sd_wrt_b_fun(Hessian_sd_wrt_b_buffer, AA[ass,:], bb[ass])

    #     get_grad_sd_wrt_xt!(grad_sd_wrt_xt_buffer, grad_sd_wrt_A_buffer, grad_A_wrt_xt_buffer, grad_sd_wrt_b_buffer, grad_b_wrt_xt_buffer)
    #     get_Hessian_sd_wrt_xt!(Hessian_sd_wrt_xt_buffer, grad_sd_wrt_A_buffer, grad_A_wrt_xt_buffer, grad_sd_wrt_b_buffer, grad_b_wrt_xt_buffer, Hessian_sd_wrt_A_buffer, Hessian_A_wrt_xt_buffer, Hessian_sd_wrt_b_buffer, Hessian_b_wrt_xt_buffer)

    #     Jsdlag .= [
    #                 [-λsd*Hessian_sd_wrt_xt_buffer zeros(6, 6); zeros(6, 6)  zeros(6, 6)]     [-grad_sd_wrt_xt_buffer; zeros(6)]
    #                 ]
    # end

    function fill_F!(F, θ, x0)
        # TODO obs_polys as parameters
        F .= 0.0 # clear
        @inbounds z = θ[z_s2i[:]]
        @inbounds λ_nom = θ[λ_nom_s2i[:]]
        get_Fnom!(F, z, x0, λ_nom)

        for t in 1:T
            xt_ind = z_s2i[1:n_x, t]
            @inbounds xt = z[xt_ind]
            # println("F ", t)
            for (i, Pe) in enumerate(ego_polys)
                m1 = length(Pe.b)

                get_grad_A_ego_wrt_xt[i](grad_A_ego_wrt_xt_buffer, xt[1:6])
                get_grad_b_ego_wrt_xt[i](grad_b_ego_wrt_xt_buffer, xt[1:6])
                # get_J_A_ego_wrt_xt_fun[i](J_A_ego_wrt_xt_buffer, xt[1:6])
                # get_J_b_ego_wrt_xt_fun[i](J_b_ego_wrt_xt_buffer, xt[1:6])
                for (j, Po) in enumerate(obs_polys)
                    (sorted_sds, sorted_ass, AA, bb) = get_sorted_sds_3d(i, length(Pe.b), j, length(Po.b), xt)
                    n_ass = min(length(sorted_ass), n_sd_slots)
                    k_map = compute_ass_ind_map(sorted_ass, n_ass)


                    for slot_i in 1:n_sd_slots
                        if slot_i < n_ass
                            # smart allocation for existing assignments
                            sd_rank = slot_i
                            k = k_map[sd_rank]
                        else
                            # copy last assignment for others
                            sd_rank = n_ass
                            k = slot_i
                        end
                        ass = sorted_ass[sd_rank]

                        sd_ind = sd_cons_s2i[k, j, i, t]
                        @inbounds λsd = θ[sd_ind]
                        
                        
                        get_grad_Ab_wrt_xt!(grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer, ass, grad_A_ego_wrt_xt_buffer, grad_b_ego_wrt_xt_buffer, m1)
                        # get_Hessian_Ab_wrt_xt!(Hessian_A_wrt_xt_buffer, Hessian_b_wrt_xt_buffer, ass, Hessian_A_ego_wrt_xt_buffer, Hessian_b_ego_wrt_xt_buffer, m1)

                        Ab = vcat(vec(AA[ass,:]), bb[ass])
                        get_grad_sd_wrt_Ab(grad_sd_wrt_Ab_buffer, Ab)
                        # get_Hessian_sd_wrt_Ab(Hessian_sd_wrt_Ab_buffer, Ab)

                        get_grad_sd_wrt_xt!(grad_sd_wrt_xt_buffer, grad_sd_wrt_Ab_buffer, grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer)
                        # get_Hessian_sd_wrt_xt!(Hessian_sd_wrt_xt_buffer, grad_sd_wrt_Ab_buffer, grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer, Hessian_sd_wrt_Ab_buffer, Hessian_A_wrt_xt_buffer, Hessian_b_wrt_xt_buffer)

                        sd_lag_buf .= [-λsd * grad_sd_wrt_xt_buffer; zeros(6)]

                        @inbounds F[xt_ind] += sd_lag_buf
                        @inbounds F[sd_ind] += sorted_sds[sd_rank]
                    end

                    # not updating sd slot memory in F is less robust
                    #sd_slot_mem[1:n_ass] = sorted_ass[1:n_ass]
                end
            end
        end
        nothing
    end

    # fill_J!
    Jnom_buf = zeros(length(Jnom_vals))

    function fill_J_vals!(J_vals, θ, x0)
        ### check Jacobian numerically
        #buf = zeros(n)
        #fill_F!(buf, θ, x0)
        #buf2 = zeros(n)
        #Jnum = spzeros(n, n)
        #for ni in 1:n
        #    wi = deepcopy(θ)
        #    wi[ni] += 1e-5
        #    fill_F!(buf2, wi, x0)
        #    #@infiltrate ni == 49
        #    Jnum[:, ni] = sparse((buf2 - buf) ./ 1e-5)
        #end
        ####

        # TODO obs_polys as parameters
        J_vals.nzval .= 1e-16 # clear
        @inbounds z = θ[z_s2i[:]]
        @inbounds λ_nom = θ[λ_nom_s2i[:]]

        get_Jnom_vals!(Jnom_buf, z, x0, λ_nom)
        J_vals .+= sparse(Jnom_rows, Jnom_cols, Jnom_buf, n, n)

        for t in 1:T
            xt_ind = z_s2i[1:n_x, t]
            @inbounds xt = z[xt_ind]
            # println("J ", t)
            for (i, Pe) in enumerate(ego_polys)
                m1 = length(Pe.b)
# @infiltrate
                get_grad_A_ego_wrt_xt[i](grad_A_ego_wrt_xt_buffer, xt[1:6])
                get_grad_b_ego_wrt_xt[i](grad_b_ego_wrt_xt_buffer, xt[1:6])
                get_Hessian_A_ego_wrt_xt[i](Hessian_A_ego_wrt_xt_buffer, xt[1:6])
                get_Hessian_b_ego_wrt_xt[i](Hessian_b_ego_wrt_xt_buffer, xt[1:6])
# @infiltrate
                for (j, Po) in enumerate(obs_polys)
                    (sorted_sds, sorted_ass, AA, bb) = get_sorted_sds_3d(i, length(Pe.b), j, length(Po.b), xt)

                    n_ass = min(length(sorted_ass), n_sd_slots)
                    k_map = compute_ass_ind_map(sorted_ass, n_ass)

                    for slot_i in 1:n_sd_slots
                        if slot_i < n_ass
                            # smart allocation for existing assignments
                            sd_rank = slot_i
                            k = k_map[sd_rank]
                        else
                            # copy last assignment for others
                            sd_rank = n_ass
                            k = slot_i
                        end
                        ass = sorted_ass[sd_rank]

                        sd_ind = sd_cons_s2i[k, j, i, t]
                        @inbounds λsd = θ[sd_ind]

                        get_grad_Ab_wrt_xt!(grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer, ass, grad_A_ego_wrt_xt_buffer, grad_b_ego_wrt_xt_buffer, m1)
                        get_Hessian_Ab_wrt_xt!(Hessian_A_wrt_xt_buffer, Hessian_b_wrt_xt_buffer, ass, Hessian_A_ego_wrt_xt_buffer, Hessian_b_ego_wrt_xt_buffer, m1)

                        Ab = vcat(vec(AA[ass,:]), bb[ass])
                        get_grad_sd_wrt_Ab(grad_sd_wrt_Ab_buffer, Ab)
                        get_Hessian_sd_wrt_Ab(Hessian_sd_wrt_Ab_buffer, Ab)

                        
                        get_grad_sd_wrt_xt!(grad_sd_wrt_xt_buffer, grad_sd_wrt_Ab_buffer, grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer)
                        get_Hessian_sd_wrt_xt!(Hessian_sd_wrt_xt_buffer, grad_sd_wrt_Ab_buffer, grad_A_wrt_xt_buffer, grad_b_wrt_xt_buffer, Hessian_sd_wrt_Ab_buffer, Hessian_A_wrt_xt_buffer, Hessian_b_wrt_xt_buffer)

                        Jsd_buf .= [grad_sd_wrt_xt_buffer; zeros(6)]
                        Jsdlag_buf .= [
                                    [-λsd*Hessian_sd_wrt_xt_buffer zeros(6, 6); zeros(6, 6)  zeros(6, 6)]     [-grad_sd_wrt_xt_buffer; zeros(6)]
                                    ]


                        # the same as code block above
                        # get_Jsdlag!(Jsdlag_buf, xt, λsd, i, ass, AA, bb, m1)
                        # get_Jsd!(Jsd_buf, xt, i, ass, AA, bb, m1)

                        @inbounds J_vals[xt_ind, xt_ind] += Jsdlag_buf[1:n_x, 1:n_x]
                        @inbounds J_vals[xt_ind, sd_ind] += Jsdlag_buf[1:n_x, n_x+1]
                        @inbounds J_vals[sd_ind, xt_ind] += Jsd_buf
                        # if maximum(J_vals) > 1e8
                        #     @infiltrate
                        #     # plot_xt(xt, Pe.A, Pe.b, Po.A, Po.b)
                        # end
                    end

                    # update sd slot memory
                    sd_slot_mem[1:n_ass] = sorted_ass[1:n_ass]
                end
            end
        end
        # check if Jacobian has exploded
        # @infiltrate
        # max_val, max_ind = findmax(J_vals)
        # min_val, min_ind = findmin(J_vals)
        # println("max", Tuple(max_ind), "=", round(max_val; sigdigits=3), " min", Tuple(min_ind), "=", round(min_val; sigdigits=3))
        # println("J[30,494]", J_vals[30,494])
        nothing
    end

    J_example = sparse(Jnom_rows, Jnom_cols, ones(length(Jnom_cols)), n, n)

    for t in 1:T
        xt_ind = z_s2i[1:n_x, t]
        @inbounds xt = z[xt_ind]

        for (i, Pe) in enumerate(ego_polys)
            for (j, Po) in enumerate(obs_polys)
                sd_ind = sd_cons_s2i[:, j, i, t]
                @inbounds J_example[xt_ind, xt_ind] .+= 1.0
                @inbounds J_example[xt_ind, sd_ind] .+= 1.0
                @inbounds J_example[sd_ind, xt_ind] .+= 1.0
            end
        end
    end

    # forced compilation and J example
    # TODO verify this is okay
    #@info "Forcing compilation by computing J example"
    #J_example = sparse(Jnom_rows, Jnom_cols, ones(length(Jnom_cols)), n, n)
    #fill_J_vals!(J_example, rand(length(θ)), rand(length(x0)))

    param = (;
        ego_polys,
        obs_polys,
        T,
        dt,
        R_cost,
        Q_cost,
        p1_max,
        p2_max,
        p3_max,
        u1_max,
        u2_max,
        u3_max,
        u4_max,
        u5_max,
        u6_max,
        n_sd_slots,
        z_s2i,
        dyn_cons_s2i,
        env_cons_s2i,
        sd_cons_s2i
    )

    (;
        fill_F!,
        fill_J_vals!,
        J_example,
        l,
        u,
        param
    )
end

# this needs to be updated for 3d
function visualize_nonsmooth_3d(x0, T, ego_polys, obs_polys; fig=Figure(), ax3=LScene(fig[1, 1], scenekw=(camera=cam3d!, show_axis=true)), θ=[], is_displaying=true, is_newsd=false)
    n_obs = length(obs_polys)
    n_ego = length(ego_polys)
    xxts = Dict()
    

    for i in 1:n_ego
        xx = x0[1:6]
        Aeb = shift_to_3D(ego_polys[i].A, ego_polys[i].b, xx)
        self_poly = ConvexPolygon3D(Aeb[1], Aeb[2])
        #plot!(ax, self_poly; color=:blue)
        plot_3D!(ax3, self_poly; color=:blue)

        for t in 1:T-1#5:1:T-1
            xxts[i, t] = Observable(x0[1:6])
            Aeb = @lift(shift_to_3D(ego_polys[i].A, ego_polys[i].b, $(xxts[i, t])))
            self_poly = @lift(ConvexPolygon3D($(Aeb)[1], $(Aeb)[2]))
            #plot!(ax, self_poly; color=:blue, linestyle=:dash)
            plot_3D!(ax3, self_poly; color=:blue, linestyle=:dash)
        end
        t = T
        xxts[i, t] = Observable(x0[1:6])
        Aeb = @lift(shift_to_3D(ego_polys[i].A, ego_polys[i].b, $(xxts[i, t])))
        self_poly = @lift(ConvexPolygon3D($(Aeb)[1], $(Aeb)[2]))
        #plot!(ax, self_poly; color=:blue, linewidth=3)
        plot_3D!(ax3, self_poly; color=:blue, linewidth=3)
    end

    colors = [:red for _ in 1:n_obs]
    for (P, c) in zip(obs_polys, colors)
        plot_3D!(ax3, P; color=c)
        #plot!(ax, P; color=c)
    end

    function update_fig(θ)
        for i in 1:n_ego
            for t in 1:T #5:5:T
                xxts[i, t][] = copy(θ[(t-1)*18+1:(t-1)*18+6])
            end
        end
    end

    if !isempty(θ)
        update_fig(θ)
    end

    if is_displaying
        display(fig)
    end

    (fig, update_fig, ax3)
end

function solve_nonsmooth_3d(prob, x0; θ0=nothing, is_displaying=true, sleep_duration=0.0)
    param = prob.param

    n_x = 12
    n_u = 6
    n_xu = n_x + n_u
    n = length(prob.l)

    if is_displaying
        (fig, update_fig) = visualize_nonsmooth_3d(x0, param.T, param.ego_polys, param.obs_polys)
    end

    # initialize
    if isnothing(θ0)
        θ0 = zeros(n)
        for t in 1:param.T
            θ0[(t-1)*n_xu+1:(t-1)*n_xu+n_x] = x0
        end
    end

    # F
    function F(n, θ, FF)
        FF .= 0.0
        prob.fill_F!(FF, θ, x0)

        if is_displaying
            update_fig(θ)
            if sleep_duration > 0
                sleep(sleep_duration)
            end
        end
        Cint(0)
    end

    # J
    J_vals = prob.J_example
    J_col = J_vals.colptr[1:end-1]
    J_len = diff(J_vals.colptr)
    J_row = J_vals.rowval
    nnz_total = length(J_vals.nzval)

    function J(n, nnz, θ, col, len, row, data)
        @assert nnz == nnz_total
        data .= 0.0

        prob.fill_J_vals!(J_vals, θ, x0)
        col .= J_col
        len .= J_len
        row .= J_row
        data .= J_vals.nzval
        Cint(0)
    end

    # force compilation
    buf = zeros(n)
    Jbuf = zeros(nnz_total)
    w = randn(length(θ0))
    F(n, w, buf)
    J(n, nnz_total, w, zero(J_col), zero(J_len), zero(J_row), Jbuf)

    # check Jacobian
    buf2 = zeros(n)
    Jrows, Jcols, _ = findnz(prob.J_example)
    Jnum = sparse(Jrows, Jcols, Jbuf)
    Jnum2 = spzeros(n, n)
    @info "Testing Jacobian accuracy numerically"
    @showprogress for ni in 1:n
       wi = copy(w)
       wi[ni] += 1e-8
       F(n, wi, buf2)
       Jnum2[:, ni] = sparse((buf2 - buf) ./ 1e-8)
       if norm(buf2 - buf)>1e-3
        @infiltrate
       end
    end
    @info "Jacobian error is $(norm(Jnum2-Jnum))"
    PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")
    status, θ, info = PATHSolver.solve_mcp(
        F,
        J,
        prob.l,
        prob.u,
        θ0;
        silent=true,
        nnz=nnz_total,
        jacobian_structure_constant=true,
        output_linear_model="no",
        preprocess=1,
        output_warnings="no",
        jacobian_data_contiguous=true,
        cumulative_iteration_limit=100_000,
        convergence_tolerance=5e-4
    )

    f_res = zeros(n)
    F(n, θ, f_res)

    if is_displaying
        display(fig)
    end

    @inbounds z = @view(θ[param.z_s2i[:]])

    (; status, info, θ, z, f_res)
end