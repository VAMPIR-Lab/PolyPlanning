struct UseOSQPSolver end
struct UsePATHSolver end

function solve_qp(::UseOSQPSolver;
                  P=nothing,
                  q=nothing,
                  A=nothing,
                  l=nothing,
                  u=nothing,
                  kwargs...)
    m = OSQP.Model()
    OSQP.setup!(m; P, q, A, l, u, kwargs...)
    ret = OSQP.solve!(m)
    ret
end

function solve_lmcp(::UsePATHSolver,
                    M,
                    q,
                    l,
                    u;
                    z=nothing)
    if isnothing(z)
        z0 = zero(q)
    else
        z0 = z
    end
    M = SparseMatrixCSC{Float64, Int32}(M)
    solve_status, z, info = PATHSolver.solve_mcp(M, q, l, u, z0, 
                               silent=true, 
                               convergence_tolerance=1e-8, 
                               cumulative_iteration_limit=200000,
                               major_iteration_limit=1000,
                               restart_limits=5,
                               lemke_rank_deficiency_iterations=1000)
    (; z, solve_status, info)
end

