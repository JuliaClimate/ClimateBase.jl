# Requires the A,B,C arrays created in `runtests.jl`
using GeoMakie, GLMakie

fig, ax = climplot(A)
@test ax isa GeoMakie.GeoAxis
