
##==============================================================================
## LEGACY SUPPORT FOR ZMQ IN CAESAR
##==============================================================================

export listSolvekeys

export _evalType

# not sure if and where this is still being used
function _evalType(pt::String)::Type
  @error "_evalType has been deprecated, use DFG serialization methods instead."
  try
    getfield(Main, Symbol(pt))
  catch ex
    io = IOBuffer()
    showerror(io, ex, catch_backtrace())
    err = String(take!(io))
    error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
  end
end



##==============================================================================
## Cannot delete until GraphProductOperations usage of this function updated
##==============================================================================

export findRelatedFromPotential

"""
    $(SIGNATURES)

Compute proposal belief on `vertid` through `fct` representing some constraint in factor graph.
Always full dimension variable node -- partial constraints will only influence subset of variable dimensions.
The remaining dimensions will keep pre-existing variable values.

Notes
- fulldim is true when "rank-deficient" -- TODO swap to false (or even float)
"""
function findRelatedFromPotential(dfg::AbstractDFG,
                                  fct::DFGFactor,
                                  varid::Symbol,
                                  N::Int,
                                  dbg::Bool=false;
                                  solveKey::Symbol=:default )
  #
  Base.depwarn("findRelatedFromPotential likely to be deprecated, use `lsf` or `productbelief(fg, variableSym, ...) instead`", :findRelatedFromPotential)

  # assuming it is properly initialized TODO
  ptsbw = evalFactor(dfg, fct, varid, solveKey=solveKey, N=N, dbg=dbg);
  # determine if evaluation is "dimension-deficient"

  # solvable dimension
  inferdim = getFactorSolvableDim(dfg, fct, varid)
  # zdim = getFactorDim(fct)
  # vdim = getVariableDim(DFG.getVariable(dfg, varid))

  # TODO -- better to upsample before the projection
  Ndim = size(ptsbw,1)
  Npoints = size(ptsbw,2)
  # Assume we only have large particle population sizes, thanks to addNode!
  manis = getManifolds(dfg, varid)
  # manis = getVariableType(DFG.getVariable(dfg, varid)).manifolds # older
  p = AMP.manikde!(ptsbw, manis)
  if Npoints != N # this is where we control the overall particle set size
      p = resample(p,N)
  end
  return (p, inferdim)
end


##==============================================================================
## TODO deprecated  
##==============================================================================



##==============================================================================
## Deprecate at v0.19
##==============================================================================


@deprecate setCliqueDrawColor!(w...;kw...) setCliqueDrawColor!(w...;kw...)

@deprecate evalFactor2(w...;kw...) evalFactor(w...;kw...)

# TreeBelief field softtype->variableType rename
function Base.getproperty(x::TreeBelief,f::Symbol)
  if f == :softtype
    Base.depwarn("`TreeBelief` field `softtype` is deprecated, use `variableType`", :getproperty)
    f = :variableType
  end
  getfield(x,f)
end

function Base.setproperty!(x::TreeBelief, f::Symbol, val)
  if f == :softtype
    Base.depwarn("`TreeBelief` field `softtype` is deprecated, use `variableType`", :getproperty)
    f = :variableType
  end
  return setfield!(x, f, convert(fieldtype(typeof(x), f), val))
end

function MetaBayesTree(tree::BayesTree)
  Base.depwarn("Graphs.jl Bayes Tree is deprecated, this constructor will be removed", :MetaBayesTree)
  # create graph from Graphs.jl adjacency_matrix
  mtree = MetaBayesTree(MetaDiGraph{Int, Float64}(MetaGraphs.SimpleDiGraph(Graphs.adjacency_matrix(tree.bt))), tree.btid, tree.frontals, tree.variableOrder, tree.buildTime)

  #deep copy over properties
  for v in tree.bt.vertices
    # set_prop!(mtree.bt, v.id, :label, deepcopy(v.label))
    set_prop!(mtree.bt, v.id, :clique, deepcopy(v))
  end

  return mtree

end

##==============================================================================
## Deprecate at v0.18 
##==============================================================================


# Keep these a bit longer

@deprecate wipeBuildNewTree!(dfg::AbstractDFG; kwargs...) resetBuildTree!(dfg; kwargs...)

@deprecate LinearConditional(N::Int=1) LinearRelative{N}(LinearAlgebra.I)

@deprecate LinearConditional(x...) LinearRelative(x...)

@deprecate PackedLinearConditional(x...) PackedLinearRelative(x...)

@deprecate extractdistribution(x) convert(SamplableBelief, x)


#
