
export LinearRelative, PackedLinearRelative

# TODO standardize
getDimension(::InstanceType{LinearRelative{N}}) where {N} = N

# new and simplified interface for both nonparametric and parametric
function (s::CalcFactor{<:LinearRelative})(z, x1, x2)
  # TODO convert to distance(distance(x2,x1),z) # or use dispatch on `-` -- what to do about `.-`
  # if s._sampleIdx < 5
  #   @info "LinearRelative" s._sampleIdx "$z" "$x1" "$x2" s.solvefor getLabel.(s.fullvariables)
  #   @info "in variables" pointer(getVal(s.fullvariables[s.solvefor])) getVal(s.fullvariables[s.solvefor])[1]
  # end
  return z .- (x2 .- x1)
end

function Base.convert(
  ::Type{<:MB.AbstractManifold},
  ::InstanceType{LinearRelative{N}},
) where {N}
  return Manifolds.TranslationGroup(N)
end

