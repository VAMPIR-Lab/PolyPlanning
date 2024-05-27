using Revise
using PolyPlanning

# generate polys
polys = PolyPlanning.gen_polys(3)
(P1, P2, P3) = polys

# plot a poly
fig = PolyPlanning.Figure()
ax = PolyPlanning.Axis(fig[1,1], aspect=PolyPlanning.DataAspect())
P = P1
PolyPlanning.plot!(ax, P)

# plot polys shrinking to a point near to the centroid
c = sum(P.V) / length(P.V)
S = [-0.999, -0.99, -0.9, -0.6, -0.3, 0]
fig = PolyPlanning.Figure()
ax = PolyPlanning.Axis(fig[1,1], aspect=PolyPlanning.DataAspect())
ShrunkPolys = PolyPlanning.GenerateShrinkingPolys(P1; c=c, S=S)
for ShrunkPoly in ShrunkPolys
    PolyPlanning.plot!(ax, ShrunkPoly)
end

# plot polys shrinking to a point near to the first vertex
c = PolyPlanning.WeightedAveragePoint(P.V, [100, 1, 1, 1])
fig = PolyPlanning.Figure()
ax = PolyPlanning.Axis(fig[1,1], aspect=PolyPlanning.DataAspect())
ShrunkPolys = PolyPlanning.GenerateShrinkingPolys(P1; c=c, S=S)
for ShrunkPoly in ShrunkPolys
    PolyPlanning.plot!(ax, ShrunkPoly)
end

