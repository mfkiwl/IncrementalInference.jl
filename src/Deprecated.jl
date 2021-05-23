
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
## Deprecate code below before v0.26
##==============================================================================

@deprecate ensureAllInitialized!(w...;kw...) initAll!(w...;kw...)

@deprecate getFactorMean(w...) IIF.getParametricMeasurement(w...)[1]

# """
#     $SIGNATURES

# Recover the mean (Gaussian) or estimate stochastic mean (non-Gaussian) value stored in a factor measurement.

# Related

# accumulateFactorMeans, solveBinaryFactorParameteric
# """
# function getFactorMean(fct::FunctorInferenceType)
#   fctt = typeof(fct)
#   error("no getFactorMean defined for $(fctt.name), has fields $(fieldnames(fctt))")
# end

# getFactorMean(fct::Normal) = [fct.μ]
# getFactorMean(fct::MvNormal) = fct.μ
# getFactorMean(fct::Union{<:BallTreeDensity,<:ManifoldKernelDensity}) = getKDEMean(fct)
# getFactorMean(fct::AliasingScalarSampler) = Statistics.mean(rand(fct,1000))

# getFactorMean(fct::DFGFactor) = getFactorMean(getFactorType(fct))

# getFactorMean(dfg::AbstractDFG, fctsym::Symbol) = getFactorMean(getFactor(dfg, fctsym))


# AMP.getManifolds(::T) where {T <: InferenceVariable} = getManifolds(getManifold(T))
# AMP.getManifolds(::Type{T}) where {T <: InferenceVariable} = getManifolds(getManifold(T))

getManifolds(::InstanceType{T}) where {T <: Union{InferenceVariable, AbstractFactor}} = getManifolds(getManifold(T))




##==============================================================================
## Deprecate code below before v0.25
##==============================================================================



# #TODO Consolidate with updateFromSubgraph_StateMachine
# function updateFromSubgraph_ParametricStateMachine(csmc::CliqStateMachineContainer)

#   # transfer results to main factor graph
#   frontsyms = getFrontals(csmc.cliq)
#   logCSM(csmc, "11, finishingCliq -- going for transferUpdateSubGraph! on $frontsyms")
#   transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frontsyms, updatePPE=false, solveKey=:parametric)

#   #solve finished change color
#   setCliqueDrawColor!(csmc.cliq, "lightblue")

#   logCSM(csmc, "Clique $(csmc.cliq.id): Finished", loglevel=Logging.Info)
#   return IncrementalInference.exitStateMachine

# end


# """
#     $(SIGNATURES)

# Multiply various full and partial dimension proposal densities.

# DevNotes
# - FIXME consolidate partial and full product AMP API, relates to #1010
# - TODO better consolidate with full dimension product
# - TODO -- reuse memory rather than rand here
# """
# function prodmultiplefullpartials(dens::Vector{BallTreeDensity},
#                                   partials::Dict{Any, Vector{BallTreeDensity}},
#                                   Ndims::Int,
#                                   N::Int,
#                                   manis::Tuple;
#                                   useExisting::Bool=false )
#   #
#   # calculate products over all dimensions, legacy proposals held in `dens` vector
#   pGM = AMP.manifoldProduct(dens, manis, Niter=1) |> getPoints

#   _partialProducts!(pGM, partials, manis, useExisting=useExisting)

#   return pGM
# end

# function _setCCWDecisionDimsConv!(ccwl::Union{CommonConvWrapper{F},
#                                               CommonConvWrapper{Mixture{N_,F,S,T}}} ) where {N_,F<:AbstractRelativeRoots,S,T}
#   #
#   # return nothing

#   p = Int[1:ccwl.xDim;]
#   ccwl.partialDims = SVector(Int32.(p)...)

#   # should be done with constructor only 
#   for thrid in 1:Threads.nthreads()
#     # length(ccwl.cpt[thrid].p) != ccwl.xDim ? resize!(ccwl.cpt[thrid].p, ccwl.xDim) : nothing
#     ccwl.cpt[thrid].p = p  # SVector(Int32[1:ccwl.xDim;]...)
#   end
#   nothing
# end




#
