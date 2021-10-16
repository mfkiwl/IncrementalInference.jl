# heatmap sampler (experimental)

@info "IncrementalInference.jl is loading tools related to Interpolations.jl."

using .Interpolations 

# only export on Requires.jl
export HeatmapGridDensity, PackedHeatmapGridDensity
export LevelSetGridNormal, PackedLevelSetGridNormal

export sampleHeatmap

##

getManifold(hgd::HeatmapGridDensity) = getManifold(hgd.densityFnc)
getManifold(lsg::LevelSetGridNormal) = getManifold(lsg.heatmap)

"""
    $SIGNATURES

Get the grid positions at the specified height (within the provided spreads)

DevNotes
- Recheck again that `thres` works right, as a way to discard beyond 3*sigma 
"""
function sampleHeatmap( roi::AbstractMatrix{<:Real},
                        x_grid::AbstractVector{<:Real}, 
                        y_grid::AbstractVector{<:Real},
                        # sigma_scale::Real=3,
                        thres::Real = 0  )
  #

  # mask the region of interest above the sampling threshold value
  mask = thres .<= roi

  idx2d = findall(mask)  # 2D indices
  pos = (v->[x_grid[v[1]],y_grid[v[2]]]).(idx2d)
  weights = (v->_roi[v[1],v[2]]).(idx2d)
  weights ./= sum(weights)

  # recast to the appropriate shape
  @cast kp[i,j] := pos[j][i]
  kp, weights
end


# TODO make n-dimensional, and later on-manifold
# TODO better standardize for heatmaps on manifolds w MKD
function fitKDE(support,
                weights,
                x_grid::AbstractVector{<:Real}, 
                y_grid::AbstractVector{<:Real};
                bw_factor::Real=0.7  )
  #
  # 1. set the bandwidth 
  x_spacing = Statistics.mean(diff(x_grid))
  y_spacing = Statistics.mean(diff(y_grid))
  kernel_ = bw_factor*0.5*(x_spacing + y_spacing) # 70% of the average spacing
  kernel_bw = [kernel_; kernel_]                  # same bw in x and y
  # fit KDE
  kde!(support, kernel_bw, weights)
end

# Helper function to construct HGD
function HeatmapGridDensity(data::AbstractMatrix{<:Real}, 
                            domain::Tuple{<:AbstractVector{<:Real},<:AbstractVector{<:Real}},
                            hist_callback::Union{<:Function, Nothing}=nothing,
                            bw_factor::Real=0.7,  # kde spread between domain points 
                            N::Int=10000  )
  #
  support_, weights_ = sampleHeatmap(data, domain..., 0)

  # constuct a pre-density from which to draw intermediate samples
  # TODO remove extraneous collect()
  density_ = fitKDE(collect(support_), weights_, domain...; bw_factor=bw_factor)
  pts_preIS, = sample(density_, N)
  
  @cast vec_preIS[j][i] := pts_preIS[i,j]
  
  # weight the intermediate samples according to interpolation of raw data
  hm = Interpolations.LinearInterpolation( domain, roi ) # interpolated heatmap
  d_scalar = Vector{Float64}( undef, length(vec_preIS) )
  
  # interpolate d_scalar for intermediate test points
  for (i,u) in enumerate(vec_preIS)
    if maximum(domain[1]) < abs(u[1]) || maximum(domain[2]) < abs(u[2]) 
      d_scalar[i] = 0.0
      continue
    end
    d_scalar[i] = hm(u...)
  end
  
  #
  weights = exp.(-d_scalar) # unscaled Gaussian
  weights ./= sum(weights)  # normalized
  
  # final samplable density object
  # TODO better standardize for heatmaps on manifolds
  bw = getBW(density_)[:,1]
  @cast pts[i,j] := vec_preIS[j][i]
  bel = kde!(collect(pts), bw, weights)
  density = ManifoldKernelDensity(TranslationGroup(Ndim(bel)), bel)

  # return `<:SamplableBelief` object
  HeatmapGridDensity(data, domain, hist_callback, bw_factor, density)
end


function LevelSetGridNormal(data::AbstractMatrix{<:Real}, 
                            domain::Tuple{<:AbstractVector{<:Real},<:AbstractVector{<:Real}},
                            level::Real,
                            sigma::Real;
                            sigma_scale::Real=3,
                            hist_callback::Union{<:Function, Nothing}=nothing,
                            bw_factor::Real=0.7,  # kde spread between domain points 
                            N::Int=10000  )
  #

  # select the support from raw data
  roi = data.-level
  # make Gaussian
  roi .^= 2
  roi .*= 0.5/(sigma^2)
  roi .-= sigma_scale^2
  roi .*= -1
  # truncate sigma_scale*sigma below zero
  #   h = heatmap;  z = measurement
  #   l = 1/2 (h-z/σ)^2
  #   masked_roi = 0 .< κ^2 - l
  
  hgd = HeatmapGridDensity(data, domain, hist_callback, bw_factor, N)

  LevelSetGridNormal(level, sigma, float(sigma_scale), hgd)
end




#