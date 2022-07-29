
##==============================================================================
## LEGACY, towards Sidecar
##==============================================================================

"""
Converter: Prior -> PackedPrior::Dict{String, Any}

FIXME see DFG #590 for consolidation with Serialization and Marshaling
"""
function convert(::Type{Dict{String, Any}}, prior::IncrementalInference.Prior)
    @error("Obsolete, use pack/unpack converters instead")
    z = convert(Type{Dict{String, Any}}, prior.Z)
    return Packed_Factor([z], "Prior")
end

"""
Converter: PackedPrior::Dict{String, Any} -> Prior

FIXME see DFG #590 for consolidation on Serialization and Marshaling
"""
function convert(::Type{<:Prior}, prior::Dict{String, Any})
    @error("Obsolete, use pack/unpack converters instead")
    # Genericize to any packed type next.
    z = prior["measurement"][1]
    z = convert(DFG.getTypeFromSerializationModule(z["distType"]), z)
    return Prior(z)
end


##==============================================================================
## Deprecate code below before v0.32
##==============================================================================


# """
#     $SIGNATURES
# Get `.factormetadata` for each CPT in CCW for a specific factor in `fg`. 
# """
# _getFMdThread(ccw::CommonConvWrapper, 
#               thrid::Int=Threads.threadid()) = ccw.cpt[thrid].factormetadata
# #
# _getFMdThread(fc::Union{GenericFunctionNodeData,DFGFactor}, 
#               thrid::Int=Threads.threadid()) = _getFMdThread(_getCCW(fc), thrid)
# #
# _getFMdThread(dfg::AbstractDFG,
#               lbl::Symbol,
#               thrid::Int=Threads.threadid()) = _getFMdThread(_getCCW(dfg, lbl), thrid)
# #


##==============================================================================
## Deprecate code below before v0.31
##==============================================================================

@deprecate initManual!(w...;kw...) initVariable!(w...;kw...)


