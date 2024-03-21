function f(z, T, goal_dir, R)
    cost = 0.0
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+3])
        ut = @view(z[(t-1)*6+4:(t-1)*6+6])
        #cost += ut'*R*ut - goal_dir'*xt[1:2]
        #cost += xt[1:2]'*xt[1:2]
        cost += -0.01*goal_dir'*xt[1:2] + ut'*R*ut
        #cost += 0.1*(xt[1:2]-goal_dir)'*(xt[1:2]-goal_dir) + ut'*R*ut
    end
    cost
end

function pointmass_dyn(x, u, dt)
    p1, p2, v1, v2 = x
    a1, a2 = u
    x + dt * [v1 + dt/2 * a1, v2 + dt/2 * a2, a1, a2]
end

function kinematic_bicycle_dyn(x, u, dt, L)
    p1, p2, θ, v = x
    δ, a = u
    #x + dt * [cos(θ)*v, sin(θ)*v, tan(δ)*v / L, a]
    x + dt * [cos(θ+δ/2)*(v+a/2), sin(θ+δ/2)*(v+a/2), δ, a]
end

function identity_dyn(x, u, dt)
    x + dt*u
end

function g_dyn(z, x0, T, dt, L)
    g = Num[] 
    x_prev = x0
    for t in 1:T
        xt = @view(z[(t-1)*6+1:(t-1)*6+3])
        ut = @view(z[(t-1)*6+4:(t-1)*6+6])
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
        xt = @view(z[(t-1)*6+1:(t-1)*6+3])
        ut = @view(z[(t-1)*6+4:(t-1)*6+6])
        #append!(g, [p1_max + xt[1], p1_max-xt[1], xt[2]-p2_min, xt[3]-2π, 2π-xt[3]])
        append!(g, [p1_max + xt[1], p1_max-xt[1], xt[2]-p2_min, u1_max-ut[1], ut[1]+u1_max, u2_max-ut[2], ut[2]+u2_max, u3_max-ut[3],ut[3]+u3_max, ])
    end
    g
end

function poly_from(x::AbstractArray{T}, angles, lengths) where T
    m = length(angles)
    A = zeros(T, m, 2)
    b = zeros(T, m)
    p1, p2, θ = x
    p = [p1,p2]
    for i in 1:m
        θi = θ + angles[i]
        #θi = angles[i]
        ai = [cos(θi), sin(θi)]
        bi = lengths[i] - ai'*p
        A[i,:] += ai
        b[i] += bi
    end
   
    A, b
end

function shift_to(A,b,x)
    p1,p2,θ = x
    R = [cos(θ) sin(θ);
         -sin(θ) cos(θ)]
    At = A*R'
    bt = b-At*x[1:2]
    At,bt
end

function verts_from(x::AbstractArray{T}, angles, lengths) where T
    m = length(angles)
    V = zeros(T, m, 2)

    A,b = poly_from(x, angles, lengths)
    for i in 1:m-1
        V[i,:] = -A[i:i+1,:]\b[i:i+1]
    end
    V[m,:] = -A[[m,1],:]\b[[m,1]]
    V
end

function gen_LP_data(A1::AbstractArray{T},b1,A2,b2) where T
    m1 = length(b1)
    m2 = length(b2)
    A = [[A1; A2] ones(T, m1+m2)]
    b = [b1; b2]
    q = [0, 0, 1.0]
    (A, b, q)
end


function g_col_single(xt, Ae, be, Ao, bo)
    sds = Dict()
    Aex,bex = shift_to(Ae, be, xt)
    AA, bb, qq = gen_LP_data(Aex,bex,Ao,bo)
    m1 = length(bex)
    m2 = length(bo)
    M = [zeros(Num, 3, 3) -AA';
         AA zeros(Num, m1+m2,m1+m2)]
    r = [qq; bb]
    all_active_inds = collect(1:m1+m2)
    Itr = powerset(all_active_inds) |> collect
    for active_inds in Itr
        length(active_inds) > 3 && continue
        assignment = [i ∈ active_inds for i in 1:m1+m2]
        try
            AA_active = collect(AA[active_inds,:])
            bb_active = collect(bb[active_inds])
            # TODO what if not unique primal? Need to resolve
            if length(active_inds) == 3
                zz = -AA_active\bb_active
                #zz = -(AA_active'*AA_active)\AA_active'*bb_active
            else
                zz = -AA_active'*((AA_active*AA_active')\bb_active)
            end
            sd = zz[3]
            sds[active_inds] = sd
        catch err
            if err isa LinearAlgebra.SingularException
                continue
            else
                @infiltrate
            end
        end
    end
    sds
end

function get_single_sd_ids(xt,Ae,be,Ao,bo,max_derivs)
    Aex, bex = shift_to(Ae,be,xt)
    AA, bb, qq = gen_LP_data(Aex,bex,Ao,bo)
    m1 = length(bex)
    m2 = length(bo)

    ret = solve_qp(UseOSQPSolver(); A=sparse(AA), l=-bb, q=qq, polish=true, verbose=false)
    #if ret.info.status_polish == -1
    #    @warn "not polished"
    #end
    primals = ret.x
    duals = -ret.y

    cons = AA*primals+bb
    I1 = duals .≥ 1e-2 .&& cons .< 1e-2
    I2 = duals .< 1e-2 .&& cons .< 1e-2
    I3 = duals .< 1e-2 .&& cons .≥ 1e-2
            
    sd = primals[3]
           
    all_inds = collect(1:m1+m2)
    weak_ind_options = powerset(all_inds[I2]) |> collect
    if length(weak_ind_options) > max_derivs
        @warn("More derivatives than accounted for! |I2| = $(sum(I2)). λ = $duals. cons = $cons")
        weak_ind_options = weak_ind_options[1:max_derivs]
    end
    if length(weak_ind_options) < max_derivs
        append!(weak_ind_options, [weak_ind_options[end] for _ in 1:max_derivs-length(weak_ind_options)])
    end

    assignments = map(weak_ind_options) do inds
        sort([all_inds[I1]; inds])
    end
end

function setup_quick(; T = 1,
                 dt = 1.0,
                L = 1.0,
                goal_dir = [0, -1.0],
                R = 0.01*I(3),
                p1_max = 3.0,
                p2_min = 0.0,
                u1_max = 1.0,
                u2_max = 1.0,
                u3_max = π/4,
                sides_per_poly=4,
                derivs_per_sd=4,
                N_polys=4,
                rng = MersenneTwister(420))
    
    Ae = Symbolics.@variables(Ae[1:sides_per_poly,1:2])[1] |> Symbolics.scalarize
    be = Symbolics.@variables(be[1:sides_per_poly])[1] |> Symbolics.scalarize
    Ao = Symbolics.@variables(Ao[1:sides_per_poly,1:2])[1] |> Symbolics.scalarize
    bo = Symbolics.@variables(bo[1:sides_per_poly])[1] |> Symbolics.scalarize

    Ae =  [0.0  -1.0;
           -1.0   0.0;
           0.0   1.0;
           1.0   0.0]
    be = [0.5, 0.15, 0.5,0.15]

    z = Symbolics.@variables(z[1:6*T])[1] |> Symbolics.scalarize
    xt = Symbolics.@variables(xt[1:3])[1] |> Symbolics.scalarize
    λsd = Symbolics.@variables(λsd)
    x0 = Symbolics.@variables(x0[1:3])[1] |> Symbolics.scalarize

    #sds = g_col_all(z, T, A1,b1,A2,b2,A3,b3, angles, lengths)
    sds = g_col_single(xt, Ae, be, Ao, bo)
    num_sd_cons = T*derivs_per_sd*N_polys
    
    cost = f(z, T, goal_dir, R)
    cons_dyn = g_dyn(z, x0, T, dt, L)
    cons_env = g_env(z, T, p1_max, p2_min, u1_max, u2_max, u3_max)

    cons_nom = [cons_dyn; cons_env]
    λ_nom = Symbolics.@variables(λ_nom[1:length(cons_nom)])[1] |> Symbolics.scalarize    
    λ_col = Symbolics.@variables(λ_col[1:num_sd_cons])[1] |> Symbolics.scalarize

    θ = [z; λ_nom; λ_col]

    lag = cost - cons_nom'*λ_nom
    grad_lag = Symbolics.gradient(lag, z)
    F_nom = [grad_lag; cons_nom; zeros(Num, num_sd_cons)]
    
    l = [fill(-Inf, length(grad_lag)); fill(-Inf, length(cons_dyn)); zeros(length(cons_env)); zeros(num_sd_cons)]
    u = [fill(+Inf, length(grad_lag)); fill(Inf, length(cons_dyn)); fill(Inf, length(cons_env)); fill(Inf, num_sd_cons)]
    n = length(l)

    Jnom = Symbolics.sparsejacobian(F_nom, θ)
    (Jnom_rows, Jnom_cols, Jnom_vals) = findnz(Jnom)
    get_Fnom! = Symbolics.build_function(F_nom, z, x0, λ_nom; expression=Val(false))[2]
    get_Jnom_vals = Symbolics.build_function(Jnom_vals, z, x0, λ_nom; expression=Val(false))[1]

    sd_grad_xt_lag_funcs = []
    sd_funcs = []

    get_lag = Dict() 
    get_sd = Dict()
    get_Jlag = Dict()
    get_Jsd = Dict()

    println(length(sds))
    @info "Generating symbolic solutions to sd calculation"
    @showprogress for (k, sd) in sds
        #λ1 = λ_col[(t-1)*2*N_polys + (e-1)*2+1]
        #λ2 = λ_col[(t-1)*2*N_polys + (e-1)*2+2]

        lag = Symbolics.gradient(-sd*λsd[1], xt; simplify=false)
        #get_lag[k] = Symbolics.build_function(lag, xt,Ae,be,Ao,bo,λsd; expression=Val(false))[1] 
        #get_sd[k] = Symbolics.build_function(sd, xt, Ae,be,Ao,bo,λsd; expression=Val(false))
        get_lag[k] = Symbolics.build_function(lag, xt,Ao,bo,λsd; expression=Val(false))[1] 
        get_sd[k] = Symbolics.build_function(sd, xt, Ao,bo,λsd; expression=Val(false))

        Jlag = Symbolics.jacobian(lag, [xt; λsd]; simplify=false)
        Jsd = Symbolics.jacobian([sd], xt; simplify=false)

        #get_Jlag[k] = Symbolics.build_function(Jlag,xt,Ae,be,Ao,bo,λsd; expression=Val(false))[1]
        #get_Jsd[k] = Symbolics.build_function(Jsd,xt,Ae,be,Ao,bo,λsd; expression=Val(false))[1]
        get_Jlag[k] = Symbolics.build_function(Jlag,xt,Ao,bo,λsd; expression=Val(false))[1]
        get_Jsd[k] = Symbolics.build_function(Jsd,xt,Ao,bo,λsd; expression=Val(false))[1]
    end

    offset = length(z) + length(λ_nom)
    function fill_F!(F,z,x0,Ae,be,polys,λ_nom,λ_col)
        F .= 0.0
        get_Fnom!(F,z,x0,λ_nom)
        for t in 1:T
            
            xt_inds = (t-1)*6+1:(t-1)*6+3
            xt = @view(z[xt_inds])
            for (e,P) in enumerate(polys)
                Ao = P.A
                bo = P.b
                λte_inds = (1:derivs_per_sd) .+ ((t-1)*N_polys*derivs_per_sd+(e-1)*derivs_per_sd)
                λte = λ_col[λte_inds]
                assignments = get_single_sd_ids(xt,Ae,be,Ao,bo,derivs_per_sd)
                for (ee,assignment) in enumerate(assignments)
                    if length(assignment) > 3
                        assignment = assignment[1:3]
                    end
                    #lag = get_lag[assignment](xt,Ae,be,Ao,bo,λte[ee]) 
                    #sd = get_sd[assignment](xt,Ae,be,Ao,bo,λte[ee])
                    lag = get_lag[assignment](xt,Ao,bo,λte[ee]) 
                    sd = get_sd[assignment](xt,Ao,bo,λte[ee])
                    F[xt_inds] .+= lag
                    F[λte_inds[ee] + offset] += sd
                end
            end
        end
        nothing
    end

    J_example = sparse(Jnom_rows, Jnom_cols, ones(length(Jnom_cols)),n,n)
    for t in 1:T
        xt_inds = (t-1)*6+1:(t-1)*6+3
        for e in 1:N_polys
            λte_inds = (1:derivs_per_sd) .+ ((t-1)*N_polys*derivs_per_sd+(e-1)*derivs_per_sd)
            J_example[xt_inds,xt_inds] .= 1.0
            for l in λte_inds
                J_example[xt_inds,l+offset] .= 1.0
                J_example[l+offset, xt_inds] .= 1.0
            end
        end
    end

    function get_J_both(z, x0, Ae,be, polys,λ_nom,λ_col)
        JJ = similar(J_example)
        JJ.nzval .= 1e-16
        vals = get_Jnom_vals(z, x0, λ_nom)
        JJ += sparse(Jnom_rows, Jnom_cols, vals, n, n)

        for t in 1:T
            xt_inds = (t-1)*6+1:(t-1)*6+3
            xt = z[xt_inds]
            for (e,P) in enumerate(polys)
                Ao = collect(P.A)
                bo = P.b
                λte_inds = (1:derivs_per_sd) .+ ((t-1)*N_polys*derivs_per_sd+(e-1)*derivs_per_sd)
                λte = λ_col[λte_inds]
                assignments = get_single_sd_ids(xt,Ae,be,Ao,bo,derivs_per_sd)
                for (ee, assignment) in enumerate(assignments)
                    if length(assignment) > 3
                        assignment = assignment[1:3]
                    end
                    #Jlag = get_Jlag[assignment](xt,Ae,be,Ao,bo,λte[ee])
                    #Jsd = get_Jsd[assignment](xt,Ae,be,Ao,bo,λte[ee])
                    Jlag = get_Jlag[assignment](xt,Ao,bo,λte[ee])
                    Jsd = get_Jsd[assignment](xt,Ao,bo,λte[ee])
                    JJ[xt_inds,xt_inds] += Jlag[1:3,1:3]
                    JJ[xt_inds,λte_inds[ee]+offset] += Jlag[1:3,end]
                    JJ[λte_inds[ee]+offset, xt_inds] .+= Jsd[:]
                end
            end
        end
        JJ
    end
    
    @info "Forcing compilation"
    n_z = length(z)
    n_nom = length(λ_nom)
    n_col = length(λ_col)
    xtr = randn(3)
    Aer = randn(sides_per_poly, 2)
    ber = randn(sides_per_poly)
    Aor = randn(sides_per_poly, 2)
    bor = randn(sides_per_poly)
    λsdr = randn()
    @showprogress for k in collect(keys(sds))
        ll = get_lag[k](xtr,Aor,bor, λsdr)
        ss = get_sd[k](xtr,Aor,bor, λsdr)
        J1 = get_Jlag[k](xtr,Aor,bor,λsdr)
        J2 = get_Jsd[k](xtr,Aor,bor,λsdr)
        #ll = get_lag[k](xtr,Aer,ber,Aor,bor, λsdr)
        #ss = get_sd[k](xtr,Aer,ber,Aor,bor, λsdr)
        #J1 = get_Jlag[k](xtr,Aer,ber,Aor,bor,λsdr)
        #J2 = get_Jsd[k](xtr,Aer,ber,Aor,bor,λsdr)
    end
    
    return (; fill_F!, 
            get_J_both,
            J_example,
            l, 
            u, 
            T,
            n_z=length(z), 
            n_nom=length(λ_nom), 
            n_col=length(λ_col), 
            N_polys,
            sides_per_poly,
            p1_max,
            p2_min)
end
function solve_quick(prob, x0, polys, Pe; θ0 = nothing)
    (; fill_F!, get_J_both, J_example, l, u, T, n_z, n_nom, n_col, N_polys, sides_per_poly, p1_max, p2_min) = prob

    @assert length(polys) == N_polys

    n = length(l)
    @assert n == n_z + n_nom + n_col

    fig = Figure()
    ax = Axis(fig[1,1], aspect=DataAspect())

    Ae,be = shift_to(Pe.A, Pe.b, x0)
    self_poly = ConvexPolygon2D(Ae,be)

    plot!(ax, self_poly; color=:blue)
    for P in polys
        plot!(ax, P; color=:red)
    end

    #lines!(ax, [-p1_max, p1_max],[p2_min,p2_min], color=:black)
    #lines!(ax, [-p1_max, -p1_max], [p2_min, 10], color=:black)
    #lines!(ax, [p1_max, p1_max], [p2_min, 10], color=:black)
    display(fig)

    if isnothing(θ0) 
        θ0 = zeros(n)
        for t in 1:T
            θ0[(t-1)*6+1:(t-1)*6+3] = x0[1:3]
        end
    end


    JJ = J_example 
    J_col = JJ.colptr[1:end-1]
    J_len = diff(JJ.colptr)
    J_row = JJ.rowval
    nnz_total = length(JJ.nzval)

    function F(n, θ, result)
        result .= 0.0
        #@inbounds z = @view(θ[1:n_z])
        #@inbounds λ_nom = @view(θ[n_z+1:n_z+n_nom])
        #@inbounds λ_col = @view(θ[n_z+n_nom+1:n_z+n_nom+n_col])
        @inbounds z = θ[1:n_z]
        @inbounds λ_nom = θ[n_z+1:n_z+n_nom]
        @inbounds λ_col = θ[n_z+n_nom+1:n_z+n_nom+n_col]
        fill_F!(result, z, x0, Pe.A, Pe.b, polys, λ_nom, λ_col)
        Cint(0)
    end
    function J(n, nnz, θ, col, len, row, data)
        @assert nnz == nnz_total
        data .= 0.0
        @inbounds z = θ[1:n_z]
        @inbounds λ_nom = θ[n_z+1:n_z+n_nom]
        @inbounds λ_col = θ[n_z+n_nom+1:n_z+n_nom+n_col]
        J_mat = get_J_both(z, x0, Pe.A, Pe.b, polys, λ_nom, λ_col)
        col .= J_col
        len .= J_len
        row .= J_row
        data .= J_mat.nzval
        Cint(0)
    end

    buf = zeros(n)
    buf2 = zeros(n)
    Jbuf = zeros(nnz_total)

    F(n, θ0, buf)
    J(n, nnz_total, θ0, zero(J_col), zero(J_len), zero(J_row), Jbuf)

    PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")
    status, θ, info = PATHSolver.solve_mcp(
         F,
         J,
         l,
         u,
         θ0;
         silent=false,
         nnz=nnz_total,
         jacobian_structure_constant = true,
         output_linear_model = "no",
         preprocess = 1,
         output_warnings = "no",
         jacobian_data_contiguous = true,
         cumulative_iteration_limit = 50_000,
         #convergence_tolerance=1e-8
        )

    for t in 1:T
        xt = θ[(t-1)*6+1:(t-1)*6+3]
        Ae,be = shift_to(Pe.A, Pe.b, xt)
        self_poly = ConvexPolygon2D(Ae,be)
        plot!(ax, self_poly; color=:blue, linestyle=:dash)
    end
    display(fig)

    @inbounds z = @view(θ[1:n_z])
    @inbounds λ_nom = @view(θ[n_z+1:n_z+n_nom])
    
    (; status, info, θ, z, λ_nom)
end
