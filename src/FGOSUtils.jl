# Factor Graph OS type utilities
#  IIF methods should direclty detect extended types from user import
# of convert in their namespace

import DistributedFactorGraphs: AbstractPointParametricEst


export getPPESuggestedAll, findVariablesNear, defaultFixedLagOnTree!


# export setSolvable!

manikde!(pts::AbstractArray{Float64,2}, vartype::InferenceVariable) = manikde!(pts, getManifolds(vartype))
manikde!(pts::AbstractArray{Float64,2}, vartype::Type{<:InferenceVariable}) = manikde!(pts, getManifolds(vartype))
manikde!(pts::AbstractArray{Float64,1}, vartype::Type{ContinuousScalar}) = manikde!(reshape(pts,1,:), getManifolds(vartype))

#getVariableOrder moved to DFG

"""
    $SIGNATURES

Return N=100 measurement samples for a factor in `<:AbstractDFG`.
"""
function getMeasurements(dfg::AbstractDFG, fsym::Symbol, N::Int=100)
  fnc = getFactorFunction(dfg, fsym)
  # getSample(fnc, N)
  Xi = (v->getVariable(dfg, v)).(getVariableOrder(dfg, fsym))
  freshSamples(fnc, N)
end

"""
    $SIGNATURES

Get graph node (variable or factor) dimension.
"""
getDimension(var::DFGVariable) = getSofttype(var).dims
getDimension(fct::DFGFactor) = getSolverData(fct).fnc.zDim

"""
    $SIGNATURES

Get the folder location where debug and solver information is recorded for a particular factor graph.
"""
getLogPath(dfg::AbstractDFG) = getSolverParams(dfg).logpath

"""
    $SIGNATURES

Append `str` onto factor graph log path as convenience function.
"""
joinLogPath(dfg::AbstractDFG, str::AbstractString) = joinpath(getLogPath(dfg), str)

# """
#     $SIGNATURES
#
# Set the `solvable` parameter for either a variable or factor.
# """
# function setSolvable!(dfg::AbstractDFG, sym::Symbol, solvable::Int)
#   if isVariable(dfg, sym)
#     getVariable(dfg, sym).solvable = solvable
#   elseif isFactor(dfg, sym)
#     getFactor(dfg, sym).solvable = solvable
#   end
#   return solvable
# end

"""
    $SIGNATURES

Get the ParametricPointEstimates---based on full marginal belief estimates---of a variable in the distributed factor graph.

DevNotes
- TODO update for manifold subgroups.

Related

getVariablePPE, setVariablePosteriorEstimates!, getVariablePPE!
"""
function calcVariablePPE(var::DFGVariable,
                         softt::InferenceVariable;
                         solveKey::Symbol=:default,
                         method::Type{MeanMaxPPE}=MeanMaxPPE  )::MeanMaxPPE
  #
  P = getKDE(var)
  manis = getManifolds(softt) # getManifolds(vnd)
  ops = buildHybridManifoldCallbacks(manis)
  Pme = getKDEMean(P) #, addop=ops[1], diffop=ops[2]
  Pma = getKDEMax(P, addop=ops[1], diffop=ops[2])
  suggested = zeros(getDimension(var))
  # TODO standardize after AMP3D
  @assert length(manis) == getDimension(var)
  for i in 1:length(manis)
    mani = manis[i]
    if mani == :Euclid
      suggested[i] = Pme[i]
    elseif mani == :Circular
      suggested[i] = Pma[i]
    else
      error("Unknown manifold to find PPE, $softt, $mani")
    end
  end
  MeanMaxPPE(solveKey, suggested, Pma, Pme, now())
end
# function calcVariablePPE!(retval::Vector{Float64},
#                           var::DFGVariable,
#                           softt::InferenceVariable;
#                           method::Type{MeanMaxPPE}=MeanMaxPPE )::Nothing
#   #
#   P = getKDE(var)
#   manis = getManifolds(softt) # getManifolds(vnd)
#   ops = buildHybridManifoldCallbacks(manis)
#   Pme = getKDEMean(P, addop=ops[1], diffop=ops[2])
#   Pma = getKDEMax(P, addop=ops[1], diffop=ops[2])
#   for i in 1:length(manis)
#     mani = manis[i]
#     if mani == :Euclid
#       retval[i] = Pme[i]
#     elseif mani == :Circular
#       retval[i] = Pma[i]
#     else
#       error("Unknown manifold to find PPE, $softt, $mani")
#     end
#   end
#   nothing
# end
# """
#     $SIGNATURES
#
# Get the ParametricPointEstimates---based on full marginal belief estimates---of a variable in the distributed factor graph.
# """
# function calcVariablePPE(var::DFGVariable,
#                          softt::InferenceVariable;
#                          method::Type{<:AbstractPointParametricEst}=MeanMaxPPE  )::Vector{Float64}
#   #
#   # vect = zeros(softt.dims)
#   mmppe = calcVariablePPE(MeanMaxPPE, var, softt, method=method)
#   return mmppe.suggested
# end


# calcVariablePPE!(retvec::Vector{Float64}, var::DFGVariable; method::Type{<:AbstractPointParametricEst}=MeanMaxPPE) = calcVariablePPE!(retvec, var, getSofttype(var), method=method)
calcVariablePPE(var::DFGVariable; method::Type{<:AbstractPointParametricEst}=MeanMaxPPE, solveKey::Symbol=:default) = calcVariablePPE(var, getSofttype(var), method=method, solveKey=solveKey)
function calcVariablePPE(dfg::AbstractDFG, sym::Symbol; method::Type{<:AbstractPointParametricEst}=MeanMaxPPE, solveKey::Symbol=:default )
  var = getVariable(dfg, sym)
  calcVariablePPE(var, getSofttype(var), method=method, solveKey=solveKey)
end






# not sure if and where this is still being used
function _evalType(pt::String)::Type
    try
        getfield(Main, Symbol(pt))
    catch ex
        io = IOBuffer()
        showerror(io, ex, catch_backtrace())
        err = String(take!(io))
        error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
    end
end


"""
    $SIGNATURES

Return interger index of desired variable element.

Example
-------
```julia
pp = RoME.Point2()
getIdx(pp, :posY) # = 2
```

Internal Notes
--------------
- uses number i < 100 for index number, and
- uses +100 offsets to track the minibatch number of the requested dimension
"""
function getIdx(pp::T, sym::Symbol, i::Int=0)::Tuple{Int, Int} where {T <: Tuple}
  # i > 99 ? error("stop") : nothing
  i-=100
  for p in pp
    i,j = getIdx(p, sym, i)
    if i > 0
      return i, j
    end
  end
  return i,-1
end
getIdx(pp::Symbol, sym::Symbol, i::Int=0)::Tuple{Int, Int} = pp==sym ? (abs(i)%100+1, div(abs(i)-100,100)) : (i-1, div(abs(i)-100,100))
function getIdx(pp::V, sym::Symbol, i::Int=0)::Tuple{Int, Int} where {V <: InferenceVariable}
  return getIdx(pp.dimtype, sym)
end



"""
    $SIGNATURES

Return `::Bool` on whether this variable has been marginalized.
"""
isMarginalized(vert::DFGVariable) = getSolverData(vert).ismargin
isMarginalized(dfg::AbstractDFG, sym::Symbol) = isMarginalized(DFG.getVariable(dfg, sym))

function setThreadModel!(fgl::AbstractDFG;
                         model=IncrementalInference.SingleThreaded)
  #
  for (key, id) in fgl.fIDs
    getSolverData(getFactor(fgl, key)).fnc.threadmodel = model
  end
  nothing
end

"""
    $SIGNATURES

Return bool on whether a certain factor has user defined multihypothesis.

Related

getMultihypoDistribution
"""
isMultihypo(fct::DFGFactor) = isa(getSolverData(fct).fnc.hypotheses, Distribution)

"""
    $SIGNATURES

Return the categorical distributed used for multihypothesis selection in a factor.

Related

isMultihypo
"""
getMultihypoDistribution(fct::DFGFactor) = getSolverData(fct).fnc.hypotheses

"""
    $SIGNATURES

Free all variables from marginalization.
"""
function dontMarginalizeVariablesAll!(fgl::G) where G <: AbstractDFG
  fgl.solverParams.isfixedlag = false
  fgl.solverParams.qfl = 9999999999
  fgl.solverParams.limitfixeddown = false
  for sym in ls(fgl)
    getSolverData(getVariable(fgl, sym)).ismargin = false
  end
  nothing
end

"""
    $SIGNATURES

Free all variables from marginalization.

Related

dontMarginalizeVariablesAll!
"""
function unfreezeVariablesAll!(fgl::AbstractDFG)
  dontMarginalizeVariablesAll!(fgl)
end

"""
    $SIGNATURES

Reset initialization flag on all variables in `::AbstractDFG`.

Notes
- Numerical values remain, but inference will overwrite since init flags are now `false`.
"""
function resetVariableAllInitializations!(fgl::AbstractDFG)
  vsyms = ls(fgl)
  for sym in vsyms
    setVariableInitialized!(getVariable(fgl, sym), :false)
  end
  nothing
end

"""
    $SIGNATURES

Enable defaults for fixed-lag-like operation by using smart message passing on the tree.

Notes:
- These are only default settings, and can be modified in each use case scenario.
- Default does not update downsolve through to leaves of the tree.
"""
function defaultFixedLagOnTree!(dfg::AbstractDFG,
                                len::Int=30;
                                limitfixeddown::Bool=true )
  #
  getSolverParams(dfg).isfixedlag = true
  getSolverParams(dfg).qfl = len
  getSolverParams(dfg).limitfixeddown = limitfixeddown
  getSolverParams(dfg)
end

"""
    $SIGNATURES

Return `::Tuple` with matching variable ID symbols and `Suggested` PPE values.

Related

getVariablePPE
"""
function getPPESuggestedAll(dfg::AbstractDFG,
                            regexFilter::Union{Nothing, Regex}=nothing )::Tuple{Vector{Symbol}, Matrix{Float64}}
  #
  # get values
  vsyms = listVariables(dfg, regexFilter) |> sortDFG
  slamPPE = map(x->getVariablePPE(dfg, x), vsyms)
  # sizes to convert to matrix
  rumax = zeros(Int, 2)
  for varr in slamPPE
    rumax[2] = length(varr)
    rumax[1] = maximum(rumax)
  end

  # populate with values
  XYT = zeros(length(slamPPE),rumax[1])
  for i in 1:length(slamPPE)
    XYT[i,1:length(slamPPE[i])] = slamPPE[i]
  end
  return (vsyms, XYT)
end

"""
    $SIGNATURES

Find and return a `::Tuple` of variables and distances to `loc::Vector{<:Real}`.

Related

findVariablesNearTimestamp
"""
function findVariablesNear(dfg::AbstractDFG,
                           loc::Vector{<:Real},
                           regexFilter::Union{Nothing, Regex}=nothing;
                           number::Int=3  )
  #

  xy = getPPESuggestedAll(dfg, regexFilter)
  dist = sum( (xy[2][:,1:length(loc)] .- loc').^2, dims=2) |> vec
  prm = (dist |> sortperm)[1:number]
  return (xy[1][prm], sqrt.(dist[prm]))
end


function convert(::Type{Tuple{BallTreeDensity,Float64}},
                 p::TreeBelief )
  @show size(p.val), size(p.bw), p.manifolds
  (AMP.manikde!(p.val, p.bw[:,1], p.manifolds), p.inferdim)
end

function convert(::Type{Tuple{BallTreeDensity,Float64}},
                 p::EasyMessage )
  @warn "EasyMessage is being deprecated, use TreeBelief instead"
  (AMP.manikde!(p.val, p.bw, p.manifolds), p.inferdim)
end

function convert(::Type{EasyMessage},
                 bel::Tuple{BallTreeDensity,Float64},
                 manifolds::T) where {T <: Tuple}
  EasyMessage(getPoints(bel[1]), getBW(bel[1])[:,1], manifolds, bel[2])
end

function convert(::Type{TreeBelief},
                 bel::Tuple{BallTreeDensity,Float64},
                 manifolds::T) where {T <: Tuple}
  @error "Dont use this convert(::Type{TreeBelief}, bel::Tuple{BallTreeDensity,Float64}, manifolds)"
  TreeBelief(getPoints(bel[1]), getBW(bel[1])[:,1:1], bel[2], ContinuousScalar(), manifolds)
end




#
