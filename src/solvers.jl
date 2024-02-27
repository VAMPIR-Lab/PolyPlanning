struct OSQPSolver end
struct PATHSolver end

function solve_qp(::OSQPSolver;
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
