# P is polygon with type ConvexPolygon2D, c can be any point inside the polygon, s is shrink/extend factor
function ConvexPolygon2DPointShrunk(P; c=0, s=0)
    A = P.A
    b = P.b
    V = P.V
    if c==0
        c = sum(V) / length(V)
    end
    b = b + s * (A * c + b)
    PolyPlanning.ConvexPolygon2D(A, b)
end

# generate polygons given a series of factor s; when s = -1, the polygon shrinks to a point, which may cause some error when ploting
function GenerateShrinkingPolys(P; c=0, S=[-0.99, -0.9, -0.8, -0.6, -0.3, 0])
    ShrunkPolys = [ConvexPolygon2DPointShrunk(P; c, s) for s in S]
    return ShrunkPolys
end

# generate a point inside the poly with given weight
function WeightedAveragePoint(V, weight)
    sum(V .* weight)/sum(weight)
end
