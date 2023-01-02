
"""
$TYPEDEF

User factor interface method for computing the residual values of factors.

Notes
- Also see #467 on API consolidation

```julia
function (cf::CalcFactor{<:LinearRelative})(res::AbstractVector{<:Real}, z, xi, xj)
  cf.variablelist
  cf.cache
  # generic on-manifold residual function 
  return distance(z, distance(xj, xi))
end
```

DevNotes
- Follow the Github project in IIF to better consolidate CCW FMD CPT CF CFM

Related 

[`CalcFactorMahalanobis`](@ref), [`CommonConvWrapper`](@ref)
"""
struct CalcFactor{
  FT <: AbstractFactor, 
  X, 
  C, 
  VT <: Tuple, 
  M <: AbstractManifold
}
  """ the interface compliant user object functor containing the data and logic """
  factor::FT
  """ what is the sample (particle) id for which the residual is being calculated """
  _sampleIdx::Int
  """ legacy support for variable values old functor residual functions.
      TBD, this is still being used by DERelative factors. """
  _legacyParams::X
  """ allow threading for either sampling or residual calculations (workaround for thread yield issue) """
  _allowThreads::Bool
  """ user cache of arbitrary type, overload the [`preambleCache`](@ref) function. NOT YET THREADSAFE """
  cache::C

  ## TODO Consolidation WIP with FactorMetadata
  # full list of variables connected to the factor
  # TODO make sure this list is of the active hypo only
  fullvariables::VT # Vector{<:DFGVariable} # FIXME change to tuple for better type stability
  # which index is being solved for?
  solvefor::Int
  manifold::M
end

# should probably deprecate the abstract type approach?
abstract type _AbstractThreadModel end

"""
$(TYPEDEF)
"""
struct SingleThreaded <: _AbstractThreadModel end
"""
$(TYPEDEF)
"""
struct MultiThreaded <: _AbstractThreadModel end

"""
$(TYPEDEF)

Main factor memory container used during inference operations -- i.e. values specific to one complete convolution operation

Notes
- CCW does not get serialized / persisted
- At writing, the assumption is there is just one CCW per factor
- Any multithreaded design needs to happens as sub-constainers inside CCW or otherwise, to carry separate memory.
- Since #467, `CalcFactor` is the only type 'seen by user' during `getSample` or function residual calculations `(cf::CalcFactor{<:MyFactor})`, s.t. `MyFactor <: AbstractRelative___`
- There also exists a `CalcFactorMahalanobis` for parameteric computations using as much of the same mechanics as possible.
- CCW is consolidated object of other previous types, FMD CPT CF CFM.

Related 

[`CalcFactor`](@ref), [`CalcFactorMahalanobis`](@ref)
"""
mutable struct CommonConvWrapper{
  T <: AbstractFactor, 
  VT <: Tuple,
  NTP <: Tuple, 
  CT,
  AM <: AbstractManifold,
  HP <: Union{Nothing, <:Distributions.Categorical{Float64, Vector{Float64}}},
  CH <: Union{Nothing, Vector{Int}},
  MT, 
  G
} <: FactorOperationalMemory
  # Basic factor topological info
  """ Values consistent across all threads during approx convolution """
  usrfnc!::T # user factor / function
  """ Consolidation from FMD, ordered tuple of all variables connected to this factor """
  fullvariables::VT
  # shortcuts to numerical containers
  """ Numerical containers for all connected variables.  Hypo selection needs to be passed 
      to each hypothesis evaluation event on user function via CalcFactor, #1321 """
  varValsAll::NTP
  """ dummy cache value to be deep copied later for each of the CalcFactor instances """
  dummyCache::CT
  # derived config parameters for this factor
  """ Factor manifold definition for frequent use (not the variables manifolds) """
  manifold::AM
  """ Which dimensions does this factor influence.  Sensitive (mutable) to both which 'solvefor index' variable and whether the factor is partial dimension """
  partialDims::Vector{<:Integer}
  """ is this a partial constraint as defined by the existance of factor field `.partial::Tuple` """
  partial::Bool
  """ coordinate dimension size of current target variable (see .fullvariables[.varidx]), TODO remove once only use under AbstractRelativeRoots is deprecated or resolved """
  xDim::Int
  """ probability that this factor is wholly incorrect and should be ignored during solving """
  nullhypo::Float64
  """ inflationSpread particular to this factor (by how much to dispurse the belief initial values before numerical optimization is run).  Analogous to stochastic search """
  inflation::Float64
  # multihypo specific field containers for recipe of hypotheses to compute
  """ multi hypothesis settings #NOTE no need for a parameter as type is known from `parseusermultihypo` """
  hypotheses::HP
  """ categorical to select which hypothesis is being considered during convolution operation """
  certainhypo::CH
  """ subsection indices to select which params should be used for this hypothesis evaluation """
  activehypo::Vector{Int}
  # buffers and indices to point numerical computations to specific memory locations
  """ user defined measurement values for each approxConv operation
      FIXME make type stable, JT should now be type stable if rest works.
      SUPER IMPORTANT, if prior=>point or relative=>tangent, see #1661 
      can be a Vector{<:Tuple} or more direct Vector{<: pointortangenttype} """
  measurement::Vector{MT}
  """ which index is being solved for in params? """
  varidx::Int
  """ Consolidation from CPT, the actual particle being solved at this moment """
  particleidx::Int
  """ working memory to store residual for optimization routines """
  res::Vector{Float64}
  """ experimental feature to embed gradient calcs with ccw """
  _gradients::G
end


#
