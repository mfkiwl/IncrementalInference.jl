
## Distributions to JSON/Packed types

DFG.packDistribution(dtr::Categorical) = PackedCategorical(; p = dtr.p)
DFG.packDistribution(dtr::Uniform) = PackedUniform(; a = dtr.a, b = dtr.b)
DFG.packDistribution(dtr::Normal) = PackedNormal(; mu = dtr.μ, sigma = dtr.σ)
DFG.packDistribution(dtr::ZeroMeanDiagNormal) = PackedZeroMeanDiagNormal(; diag = dtr.Σ.diag)
DFG.packDistribution(dtr::ZeroMeanFullNormal) = PackedZeroMeanFullNormal(; cov = dtr.Σ.mat[:])
DFG.packDistribution(dtr::DiagNormal) = PackedDiagNormal(; mu = dtr.μ, diag = dtr.Σ.diag)
DFG.packDistribution(dtr::FullNormal) = PackedFullNormal(; mu = dtr.μ, cov = dtr.Σ.mat[:])
DFG.packDistribution(dtr::Rayleigh) = PackedRayleigh(; sigma = dtr.σ)

## Unpack JSON/Packed to Distribution types

DFG.unpackDistribution(dtr::PackedCategorical) = Categorical(dtr.p ./ sum(dtr.p))
DFG.unpackDistribution(dtr::PackedUniform) = Uniform(dtr.a, dtr.b)
DFG.unpackDistribution(dtr::PackedNormal) = Normal(dtr.mu, dtr.sigma)
function DFG.unpackDistribution(dtr::PackedZeroMeanDiagNormal)
  return MvNormal(LinearAlgebra.Diagonal(map(abs2, sqrt.(dtr.diag))))
end # sqrt.(dtr.diag)
function DFG.unpackDistribution(dtr::PackedZeroMeanFullNormal)
  d = round(Int, sqrt(size(dtr.cov)[1]))
  return MvNormal(reshape(dtr.cov, d, d))
end
DFG.unpackDistribution(dtr::PackedDiagNormal) = MvNormal(dtr.mu, sqrt.(dtr.diag))
function DFG.unpackDistribution(dtr::PackedFullNormal)
  return MvNormal(dtr.mu, reshape(dtr.cov, length(dtr.mu), :))
end
DFG.unpackDistribution(dtr::PackedRayleigh) = Rayleigh(dtr.sigma)

