
export EuclidDistance, PackedEuclidDistance



"""
$(TYPEDEF)

Default linear offset between two scalar variables.
"""
struct EuclidDistance{T <: SamplableBelief} <: AbstractRelativeMinimize
  Z::T
end

EuclidDistance(::UniformScaling=LinearAlgebra.I) = EuclidDistance(Normal())

# consider a different group?
getManifold(::InstanceType{EuclidDistance}) = TranslationGroup(1)
getDimension(::InstanceType{<:EuclidDistance}) = 1


# new and simplified interface for both nonparametric and parametric
(s::CalcFactor{<:EuclidDistance})(z, x1, x2) = z .- norm(x2 .- x1)


Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{EuclidDistance}) = Manifolds.TranslationGroup(1)


"""
$(TYPEDEF)
Serialization type for `EuclidDistance` binary factor.
"""
mutable struct PackedEuclidDistance <: AbstractPackedFactor
  _type::String
  Z::String
end

function convert(::Type{PackedEuclidDistance}, d::EuclidDistance)
  PackedEuclidDistance( "/application/JuliaLang/PackedSamplableBelief",
                        convert(String, d.Z) ) # TODO, PackedSamplableBelief
end

function convert(::Type{<:EuclidDistance}, d::PackedEuclidDistance)
  EuclidDistance(convert(SamplableBelief, d.Z))
end







#