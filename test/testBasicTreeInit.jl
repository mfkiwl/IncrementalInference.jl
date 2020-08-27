
using Test
using IncrementalInference


@testset "basic test for tree initialization functionality" begin

# small canonical factor graph, without graphinit
fg = generateCanonicalFG_CaesarRing1D(graphinit=false)
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true

# temporary
# getSolverParams(fg).dbg = true

# do all init on tree as part of solve
# getSolverParams(fg).drawtree = true
tree, smt, hist = solveTree!(fg, verbose=true)


end


#
