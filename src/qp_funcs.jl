using Symbolics
using SparseArrays
using OSQP

function solve_qp_2dirs(Q,q,A,l,u)
    model = OSQP.Model()
    OSQP.setup!(model; P=sparse(Q),q,A=sparse(A),l,u,polish=true,verbose=false)
    ret = OSQP.solve!(model)
    [ret.info.obj_val, ret.info.obj_val]
end

@register_array_symbolic solve_qp_2dirs(Q::AbstractMatrix, q::AbstractVector, A::AbstractMatrix, l::AbstractVector, u::AbstractVector) begin
    size=(2,)
    eltype=eltype(Q)
end


