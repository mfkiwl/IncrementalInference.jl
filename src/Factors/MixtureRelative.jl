
export MixtureRelative, PackedMixtureRelative


"""
$(TYPEDEF)

Define a categorical mixture of (relative) likelihood beliefs between any two variables.
"""
struct MixtureRelative{N, F<:FunctorInferenceType, S, T<:Tuple} <: AbstractRelative
  mechanics::F
  components::NamedTuple{S,T}
  diversity::Distributions.Categorical
  dims::Int
  labels::Vector{Int}
end


# NamedTuple{_defaultNamesMixtures(N)}((z...,)), c

MixtureRelative(f::Type{F},
                z::NamedTuple{S,T}, 
                c::Distributions.DiscreteNonParametric ) where {F<:FunctorInferenceType, S, T} = MixtureRelative{length(z),F,S,T}(f(LinearAlgebra.I), z, c, size( rand(z[1],1), 1), zeros(Int, 0))
#
MixtureRelative(f::F,
                z::NamedTuple{S,T}, 
                c::Distributions.DiscreteNonParametric ) where {F<:FunctorInferenceType, S, T} = MixtureRelative{length(z),F,S,T}(f, z, c, size( rand(z[1],1), 1), zeros(Int, 0))
#
MixtureRelative(f::Union{F,Type{F}},z::NamedTuple{S,T}, c::AbstractVector{<:Real}) where {F<:FunctorInferenceType,S,T} = MixtureRelative(f, z, Categorical([c...]) )
MixtureRelative(f::Union{F,Type{F}},z::NamedTuple{S,T}, c::NTuple{N,<:Real}) where {N,F<:FunctorInferenceType,S,T} = MixtureRelative(f, z, [c...] )
MixtureRelative(f::Union{F,Type{F}},z::AbstractVector{<:SamplableBelief}, c::Union{<:Distributions.DiscreteNonParametric, <:AbstractVector{<:Real}, <:NTuple{N,<:Real}} ) where {F <: FunctorInferenceType, N} = MixtureRelative(f,NamedTuple{_defaultNamesMixtures(length(z))}((z...,)), c)
MixtureRelative(f::Union{F,Type{F}},z::Tuple, c::Union{<:Distributions.DiscreteNonParametric, <:AbstractVector{<:Real}, <:NTuple{N,<:Real}} ) where {F<:FunctorInferenceType, N} = MixtureRelative(f,NamedTuple{_defaultNamesMixtures(length(z))}(z), c)


function Base.resize!(mp::Union{MixturePrior, MixtureRelative}, s::Int)
  resize!(mp.labels, s)
end

# TODO make in-place memory version
function getSample(s::Union{MixturePrior, MixtureRelative}, N::Int=1)
  #out memory should be right size first
  (length(s.labels) != N) && resize!(s, N)
  s.labels .= rand(s.diversity, N)
  smpls = Array{Float64,2}(undef,s.dims,N)
  for i in 1:N
    mixComponent = s.components[s.labels[i]]
    smpls[:,i] = rand(mixComponent,1)
  end
  (smpls, s.labels)
end



(s::MixtureRelative)( res::AbstractArray{<:Real},
                      userdata::FactorMetadata,
                      idx::Int,
                      meas::Tuple,
                      X... ) = s.mechanics(res, userdata, idx, meas, X...)
#






"""
$(TYPEDEF)

Serialization type for `MixtureRelative`.
"""
mutable struct PackedMixtureRelative <: PackedInferenceType
  N::Int
  # store the packed type for later unpacking
  F_::String 
  S::Vector{String}
  components::Vector{String}
  diversity::String
end

function convert(::Type{PackedMixtureRelative}, obj::MixtureRelative{N,F,S,T}) where {N,F,S,T}
  allcomp = String[]
  for val in obj.components
    push!(allcomp, string(val))
  end
  # pm = DFG.convertPackedType(obj.mechanics)
  pm = convert(DFG.convertPackedType(obj.mechanics), obj.mechanics)
  sT = string(typeof(pm))
  PackedMixtureRelative(N, sT, string.(collect(S)), allcomp, string(obj.diversity))
end
function convert(::Type{MixtureRelative}, obj::PackedMixtureRelative)
  N = obj.N
  F1 = getfield(Main, Symbol(obj.F_))
  S = (Symbol.(obj.S)...,)
  F2 = DFG.convertStructType(F1) # typeof(obj.mechanics))
  # T = (extractdistribution.(obj.components)...,) |> typeof
  # mechanics = F2(LinearAlgebra.I)
  components = extractdistribution.(obj.components)
  diversity = extractdistribution(obj.diversity)
  # MixtureRelative{N,F2,S,T}(mechanics, components, diversity)
  tupcomp = (components...,)
  ntup = NamedTuple{S,typeof(tupcomp)}(tupcomp)
  MixtureRelative(F2, ntup, diversity)
end








  #
