
import DistributedFactorGraphs: getVariableType

"""
    CliqStatus
Clique status message enumerated type with status.
"""
@enum CliqStatus NULL NO_INIT INITIALIZED UPSOLVED MARGINALIZED DOWNSOLVED UPRECYCLED ERROR_STATUS


# Used for UPWARD_DIFFERENTIAL, UPWARD_COMMON, DOWNWARD_COMMON marginalized types
abstract type MessagePassDirection end
struct UpwardPass <: MessagePassDirection end
struct DownwardPass <: MessagePassDirection end

abstract type MessageType end
struct NonparametricMessage <: MessageType end
struct ParametricMessage <: MessageType end


const SamplableBelief = Union{Distributions.Distribution, KernelDensityEstimate.BallTreeDensity, AliasingScalarSampler, FluxModelsDistribution}

abstract type PackedSamplableBelief end

#Supported types for parametric
const ParametricTypes = Union{Normal, MvNormal}


"""
    $TYPEDEF

INTERMEDIATE DATA STRUCTURE DURING REFACTORING.

Representation of the belief of a single variable.

Notes:
- we want to send the joint, this is just to resolve consolidation #459 first.
- Long term objective is single joint definition, likely called `LikelihoodMessage`.
"""
struct TreeBelief{T <: InferenceVariable}
  val::Array{Float64,2}
  bw::Array{Float64,2}
  inferdim::Float64
  # see DFG #603, variableType defines the domain and manifold as well as group operations for a variable in the factor graph
  variableType::T
  # TODO -- DEPRECATE
  manifolds::Tuple{Vararg{Symbol}} # NOTE added during #459 effort
  # only populated during up as solvableDims for each variable in clique, #910
  solvableDim::Float64 
end
TreeBelief( p::BallTreeDensity,
            inferdim::Real=0,
            variableType::T=ContinuousScalar(),
            manifolds=getManifolds(variableType),
            solvableDim::Real=0) where {T <: InferenceVariable} = TreeBelief{T}(getPoints(p), getBW(p), inferdim, variableType, manifolds, solvableDim)

TreeBelief( val::Array{Float64,2},
            bw::Array{Float64,2},
            inferdim::Real=0,
            variableType::T=ContinuousScalar(),
            manifolds=getManifolds(variableType),
            solvableDim::Real=0) where {T <: InferenceVariable} = TreeBelief{T}(val, bw, inferdim, variableType, manifolds, solvableDim)

function TreeBelief(vnd::VariableNodeData, solvDim::Real=0)
  TreeBelief( vnd.val, vnd.bw, vnd.inferdim, getVariableType(vnd), getManifolds(vnd), solvDim )
end

TreeBelief(vari::DFGVariable, solveKey::Symbol=:default; solvableDim::Real=0) = TreeBelief( getSolverData(vari, solveKey) , solvableDim)

getVariableType(tb::TreeBelief) = tb.variableType

getManifolds(treeb::TreeBelief) = getManifolds(treeb.variableType)

function compare(t1::TreeBelief, t2::TreeBelief)
  TP = true
  TP = TP && norm(t1.val - t2.val) < 1e-5
  TP = TP && norm(t1.bw - t2.bw) < 1e-5
  TP = TP && abs(t1.inferdim - t2.inferdim) < 1e-5
  TP = TP && t1.variableType == t2.variableType
  TP = TP && abs(t1.solvableDim - t2.solvableDim) < 1e-5
  return TP
end

"""
  $(TYPEDEF)
Belief message for message passing on the tree.  This should be considered an incomplete joint probility.

Notes:
- belief -> Dictionary of [`TreeBelief`](@ref)
- variableOrder -> Ordered variable id list of the seperators in cliqueLikelihood
- cliqueLikelihood -> marginal distribution (<: `SamplableBelief`) over clique seperators.
- Older names include: productFactor, Fnew, MsgPrior, LikelihoodMessage

DevNotes:
- Used by both nonparametric and parametric.
- Objective for parametric case: `MvNormal(μ=[:x0;:x2;:l5], Σ=[+ * *; * + *; * * +])`.
- Part of the consolidation effort, see #459.
- Better conditioning for joint structure in the works using deconvolution, see #579, #635.
  - TODO confirm why <: Singleton.

  $(TYPEDFIELDS)
"""
mutable struct LikelihoodMessage{T <: MessageType} <: AbstractPrior
  status::CliqStatus
  belief::Dict{Symbol, TreeBelief}
  variableOrder::Vector{Symbol}
  cliqueLikelihood::Union{Nothing,SamplableBelief}
  msgType::T
  hasPriors::Bool
  # this is different from belief[].inferdim, as the total available infer dims remaining during down msgs -- see #910
  childSolvDims::Dict{Int, Float64} 
end


LikelihoodMessage(; status::CliqStatus=NULL,
                    beliefDict::Dict=Dict{Symbol, TreeBelief}(),
                    variableOrder::Vector{Symbol}=Symbol[],
                    cliqueLikelihood=nothing,
                    msgType::T=NonparametricMessage(),
                    hasPriors::Bool=true,
                    childSolvDims::Dict{Int, Float64}=Dict{Int, Float64}(), 
                    ) where {T <: MessageType} =
        LikelihoodMessage{T}(status, beliefDict, variableOrder, cliqueLikelihood, msgType, hasPriors, childSolvDims)
#

function Base.show(io::IO, ::MIME"text/plain", msg::LikelihoodMessage)
  t = typeof(msg)
  fields = fieldnames(t)
  nf = nfields(msg)

  for f in fields
    printstyled(io, f,": ", color=:blue)
    show(io, getproperty(msg, f))
    println(io)
  end

end

function compare( l1::LikelihoodMessage,
                  l2::LikelihoodMessage;
                  skip::Vector{Symbol}=[] )
  #
  TP = true
  TP = TP && l1.status == l2.status
  TP = TP && l1.variableOrder == l2.variableOrder
  TP = TP && l1.msgType == l2.msgType
  TP = TP && l1.cliqueLikelihood |> typeof == l2.cliqueLikelihood |> typeof
  for (k,v) in l1.belief
    TP = TP && haskey(l2.belief, k)
    TP = TP && compare(v, l2.belief[k])
  end
end

==(l1::LikelihoodMessage,l2::LikelihoodMessage) = compare(l1,l2)



## =========================================================================================
## DEPRECATE BELOW AS ABLE
## =========================================================================================


# figure out how to deprecate (not critical at the moment)
# used in RoMEPlotting
const TempUpMsgPlotting = Dict{Symbol,Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}}


"""
$(TYPEDEF)

TODO TO BE DEPRECATED
"""
mutable struct PotProd
    Xi::Symbol # Int
    prev::Array{Float64,2}
    product::Array{Float64,2}
    potentials::Array{BallTreeDensity,1}
    potentialfac::Vector{Symbol}
end

"""
$(TYPEDEF)

TODO TO BE DEPRECATED
"""
mutable struct CliqGibbsMC
    prods::Array{PotProd,1}
    lbls::Vector{Symbol}
    CliqGibbsMC() = new()
    CliqGibbsMC(a,b) = new(a,b)
end

"""
$(TYPEDEF)

TODO TO BE DEPRECATED
"""
mutable struct DebugCliqMCMC
  mcmc::Union{Nothing, Array{CliqGibbsMC,1}}
  outmsg::LikelihoodMessage
  outmsglbls::Dict{Symbol, Symbol} # Int
  priorprods::Vector{CliqGibbsMC}
  DebugCliqMCMC() = new()
  DebugCliqMCMC(a,b,c,d) = new(a,b,c,d)
end





#
