# init utils for tree based inference

export resetCliqSolve!
export addLikelihoodsDifferential!


## =============================================================================
# short preamble funcions
## =============================================================================



convert(::Type{BallTreeDensity}, src::TreeBelief) = manikde!(src.val, src.bw[:,1], src.softtype)

"""
    $(SIGNATURES)

Construct a BallTreeDensity KDE object from an IIF.TreeBelief object.

Related

manikde!, getKDE, getKDEMax, getKDEMean, TreeBelief
"""
function kde!(em::TreeBelief)
  # return AMP.manikde!(em.val, em.bw, em.manifolds)
  return convert(BallTreeDensity, em)
end
manikde!(em::TreeBelief) = convert(BallTreeDensity, em)




## =============================================================================
# helper functions for tree message channels
## =============================================================================


"""
    $SIGNATURES

Reset the state of all variables in a clique to not initialized.

Notes
- resets numberical values to zeros.

Dev Notes
- TODO not all kde manifolds will initialize to zero.
- FIXME channels need to be consolidated
"""
function resetCliqSolve!( dfg::AbstractDFG,
                          treel::AbstractBayesTree,
                          cliq::TreeClique;
                          solveKey::Symbol=:default)
  #
  cda = getCliqueData(cliq)
  vars = getCliqVarIdsAll(cliq)
  for varis in vars
    resetVariable!(dfg, varis, solveKey=solveKey)
  end
  # TODO remove once consolidation with upMsgs is done
  putCliqueMsgUp!(cda, LikelihoodMessage() )

  # cda.dwnMsg = LikelihoodMessage()
  putCliqueInitMsgDown!(cda, LikelihoodMessage() )

  setCliqueStatus!(cliq, :null)
  setCliqDrawColor(cliq, "")
  return nothing
end

function resetCliqSolve!( dfg::AbstractDFG,
                          treel::AbstractBayesTree,
                          frt::Symbol;
                          solveKey::Symbol=:default  )
  #
  resetCliqSolve!(dfg, treel, getClique(treel, frt), solveKey=solveKey)
end



## =============================================================================
# helper functions to add tree messages to subgraphs
## =============================================================================



function updateSubFgFromDownMsgs!(sfg::G,
                                  dwnmsgs::LikelihoodMessage,
                                  seps::Vector{Symbol} ) where G <: AbstractDFG
  #
  # sanity check basic Bayes (Junction) tree property
  # length(setdiff(keys(dwnmsgs), seps)) == 0 ? nothing : error("updateSubFgFromDownMsgs! -- separators and dwnmsgs not consistent")

  # update specific variables in sfg from msgs
  for (key,beldim) in dwnmsgs.belief
    if key in seps
      setValKDE!(sfg, key, manikde!(beldim.val,beldim.bw[:,1],getManifolds(beldim.softtype)), false, beldim.inferdim)
    end
  end

  nothing
end



function generateMsgPrior(belief_::TreeBelief, ::NonparametricMessage)
  kdePr = manikde!(belief_.val, belief_.bw[:,1], getManifolds(belief_.softtype))
  MsgPrior(kdePr, belief_.inferdim)
end

function generateMsgPrior(belief_::TreeBelief, ::ParametricMessage)
  msgPrior = if size(belief_.val, 2) == 1 && size(belief_.val, 1) == 1
    MsgPrior(Normal(belief_.val[1], sqrt(belief_.bw[1])), belief_.inferdim)
  elseif size(belief_.val, 2) == 1 && 1 != size(belief_.val, 1)
    mvnorm = createMvNormal(belief_.val[:,1], belief_.bw)
    mvnorm !== nothing ? nothing : (return DFGFactor[])
    MsgPrior(mvnorm, belief_.inferdim)
  end
  return msgPrior
end


"""
    $SIGNATURES

Build from a `LikelihoodMessage` a temporary distributed factor graph object containing differential
information likelihood factors based on values in the messages.

Notes
- Modifies tfg argument by adding `:UPWARD_DIFFERENTIAL` factors.

DevNotes
- Initial version which only works for Pose2 and Point2 at this stage.
"""
function addLikelihoodsDifferential!( msgs::LikelihoodMessage, 
                                      cliqSubFG::AbstractDFG,
                                      tfg::AbstractDFG=initfg() )
  # create new local dfg and add all the variables with data
  listVarByDim = Symbol[]
  for (label, val) in msgs.belief
    push!(listVarByDim, label)
    if !exists(tfg, label)
      addVariable!(tfg, label, val.softtype)
      @debug "New variable added to subfg" _group=:check_addLHDiff #TODO JT remove debug. 
    end
    initManual!(tfg, label, manikde!(val))
  end

  # list all variables in order of dimension size
  alreadylist = Symbol[]
  listDims = getDimension.(getVariable.(tfg,listVarByDim))
  per = sortperm(listDims, rev=true)
  listVarDec = listVarByDim[per]
  listVarAcc = reverse(listVarDec)
  # add all differential factors (without deconvolution values)
  for sym1_ in listVarDec
    push!(alreadylist, sym1_)
    for sym2_ in setdiff(listVarAcc, alreadylist)
      nfactype = selectFactorType(tfg, sym1_, sym2_)
      # assume default helper function # buildFactorDefault(nfactype)
      nfct = nfactype()
      afc = addFactor!(tfg, [sym1_;sym2_], nfct, graphinit=false, tags=[:DUMMY;])
      # calculate the general deconvolution between variables
      pts = solveFactorMeasurements(tfg, afc.label)
      newBel = manikde!(pts[1], getManifolds(nfactype))
      # replace dummy factor with real deconv factor using manikde approx belief measurement
      fullFct = nfactype(newBel)
      deleteFactor!(tfg, afc.label)
      addFactor!( cliqSubFG, [sym1_;sym2_], fullFct, graphinit=false, tags=[:LIKELIHOODMESSAGE; :UPWARD_DIFFERENTIAL] )
    end
  end

  return tfg
end
# default verbNoun API spec (dest, src)
addLikelihoodsDifferential!(subfg::AbstractDFG, msgs::LikelihoodMessage) = addLikelihoodsDifferential!(msgs, subfg)

"""
    $SIGNATURES
Place a single message likelihood prior on the highest dimension variable with highest connectivity in existing subfg.
"""
function addLikelihoodPriorCommon!( subfg::AbstractDFG, 
                                    msgs::LikelihoodMessage;
                                    tags::Vector{Symbol}=Symbol[])
  #
  # find max dimension variable, which also has highest biadjacency

  len = length(msgs.belief)
  dims = Vector{Int}(undef, len)
  syms = Vector{Symbol}(undef, len)
  biAdj = Vector{Int}(undef, len)
  # TODO, not considering existing priors for MsgPrior placement at this time
  # priors = Vector{Int}(undef, len)
  i = 0
  for (label, val) in msgs.belief
    i += 1
    dims[i] = getDimension(val.softtype)
    syms[i] = label
    biAdj[i] = ls(subfg, label) |> length
  end
  # work only with highest dimension variable
  maxDim = maximum(dims)
  dimMask = dims .== maxDim
  mdAdj = biAdj[dimMask]
  pe = sortperm(mdAdj, rev=true) # decending
  topCandidate = (syms[dimMask])[pe][1]

  # get prior for top candidate
  msgPrior = generateMsgPrior(msgs.belief[topCandidate], msgs.msgType)

  # get ready
  tags__ = union(Symbol[:LIKELIHOODMESSAGE;:UPWARD_COMMON], tags)
  # finally add the single AbstractPrior from LikelihoodMessage
  addFactor!(subfg, [topCandidate], msgPrior, graphinit=false, tags=tags__)
end


"""
    $SIGNATURES

Special function to add a few variables and factors to the clique sub graph required for downward solve in CSM.

Dev Notes
- There is still some disparity on whether up and down solves of tree should use exactly the same subgraph...  'between for up and frontal connected for down'
"""
function addDownVariableFactors!( dfg::AbstractDFG,
                                  subfg::InMemoryDFGTypes,
                                  cliq::TreeClique,
                                  logger=ConsoleLogger();
                                  solvable::Int=1  )
  #
  # determine which variables and factors needs to be added
  currsyms = ls(subfg)
  allclsyms = getCliqVarsWithFrontalNeighbors(dfg, cliq, solvable=solvable)
  newsyms = setdiff(allclsyms, currsyms)
  with_logger(logger) do
    @info "addDownVariableFactors!, cliq=$(cliq.index), newsyms=$newsyms"
  end
  frtls = getCliqFrontalVarIds(cliq)
  with_logger(logger) do
    @info "addDownVariableFactors!, cliq=$(cliq.index), frtls=$frtls"
  end
  allnewfcts = union(map(x->findFactorsBetweenFrom(dfg,union(currsyms, newsyms),x), frtls)...)
  newfcts = setdiff(allnewfcts, lsf(subfg))
  with_logger(logger) do
    @info "addDownVariableFactors!, cliq=$(cliq.index), newfcts=$newfcts, allnewfcts=$allnewfcts"
  end

  #TODO solvable?
  DFG.mergeGraph!(subfg, dfg, newsyms, newfcts)

  return newsyms, newfcts
end


"""
    $SIGNATURES

Modify the `subfg::AbstractDFG` to include `msgs` as priors that are used
during clique inference.

Notes
- May be used initialization or inference, in both upward and downward directions.
- `msgs` are identified by variable label `::Symbol`, and my consist of multiple beliefs.
- Message sets from different cliques are identified by clique id `::Int`.
- assume lower limit on number of particles is 5.
- messages from children stored in vector or dict.

DevNotes
- TODO Split dispatch on `dir`, rather than internal `if` statement.

Related

`deleteMsgFactors!`
"""
function addMsgFactors!(subfg::AbstractDFG,
                        msgs::LikelihoodMessage,
                        dir::Type{<:MessagePassDirection};
                        tags::Vector{Symbol}=Symbol[])
  #
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  # TODO, expand -- this deconv approach only works for NonparametricMessage at this time.
  if getSolverParams(subfg).useMsgLikelihoods && dir == UpwardPass && msgs.msgType isa NonparametricMessage
    if 0 < msgs.belief |> length
      # currently only works for nonparametric
      addLikelihoodsDifferential!(subfg, msgs)          # :UPWARD_DIFFERENTIAL
      ## FIXME, only do this if there are priors below
      if msgs.hasPriors
        prFcts = addLikelihoodPriorCommon!(subfg, msgs)   # :UPWARD_COMMON
      end
    end
  else
    svars = DFG.listVariables(subfg)
    tags__ = union(Symbol[:LIKELIHOODMESSAGE;], tags)
    dir == DownwardPass ? push!(tags__, :DOWNWARD_COMMON) : nothing
    for (msym, belief_) in msgs.belief
      if msym in svars
        msgPrior = generateMsgPrior(belief_, msgs.msgType)
        fc = addFactor!(subfg, [msym], msgPrior, graphinit=false, tags=tags__)
        push!(msgfcts, fc)
      end
    end
  end
  return msgfcts
end

function addMsgFactors!(subfg::AbstractDFG,
                        allmsgs::Dict{Int,LikelihoodMessage},
                        dir::Type{<:MessagePassDirection};
                        tags::Vector{Symbol}=Symbol[] )
  #
  allfcts = DFGFactor[]
  for (cliqid, msgs) in allmsgs
    # do each dict in array separately
    newfcts = addMsgFactors!(subfg, msgs, dir, tags=tags)
    union!( allfcts, newfcts )
  end
  return allfcts
end


"""
    $SIGNATURES

Delete from the subgraph`::AbstractDFG` prior belief `msgs` that could/would be used
during clique inference.

DevNotes
- TODO make `::Vector{Symbol}` version.

Related

`addMsgFactors!`
"""
function deleteMsgFactors!( subfg::AbstractDFG,
                            fcts::Vector{DFGFactor} )
  #
  for fc in fcts
    deleteFactor!(subfg, fc.label)
  end
end
# deleteMsgFactors!(::LightDFG{SolverParams,DFGVariable,DFGFactor}, ::Array{DFGFactor{CommonConvWrapper{MsgPrior{BallTreeDensity}},1},1})

function deleteMsgFactors!(subfg::AbstractDFG, 
                           tags::Vector{Symbol}=[:LIKELIHOODMESSAGE])
  # remove msg factors that were added to the subfg
  facs = lsf(subfg, tags=tags)
  deleteFactor!.(subfg, facs)
  return facs
end

#TODO JT can be removed, used as sanity check
function removeSeparatorPriorsFromSubgraph!(cliqSubFg::AbstractDFG, cliq::TreeClique)
  cliqSeparatorVarIds = getCliqSeparatorVarIds(cliq)
  priorIds = Symbol[]
  for v in cliqSeparatorVarIds
    facs = getNeighbors(cliqSubFg, v)
    for f in facs
      isprior = length(getFactor(cliqSubFg, f)._variableOrderSymbols) == 1
      isprior && push!(priorIds, f)
      isprior && DFG.deleteFactor!(cliqSubFg, f)
    end
  end
  return priorIds
end



## =============================================================================
## Prepare Clique Up or Down Msgs
## =============================================================================


#TODO rename to just prepCliqueMsgUp
"""
    $SIGNATURES

Prepare the upward inference messages from clique to parent and return as `Dict{Symbol}`.

Notes
- Does not require tree message likelihood factors in subfg.
- Also see #579 regarding elimited likelihoods and priors.

DevNotes
- set `msgs.hasPriors=true` only if a prior occurred here or lower down in tree branch. 
"""
function prepCliqueMsgUpConsolidated( subfg::AbstractDFG,
                                      cliq::TreeClique,
                                      status::CliqStatus=getCliqueStatus(cliq);
                                      logger=ConsoleLogger(),
                                      duplicate::Bool=true )
  #
  # get the current clique status
  sdims = getCliqVariableMoreInitDims(subfg, cliq)

  # construct init's up msg to place in parent from initialized separator variables
  hasPriors = 0 < (lsfPriors(subfg) |> length)
  msg = LikelihoodMessage(status=status, hasPriors=hasPriors)
  seps = getCliqSeparatorVarIds(cliq)
  with_logger(logger) do
    @info "$(now()), prepCliqInitMsgsUp, seps=$seps, sdims=$sdims"
  end
  for vid in seps
    var = DFG.getVariable(subfg, vid)
    var = duplicate ? deepcopy(var) : var
    if isInitialized(var)
      msg.belief[var.label] = TreeBelief(var, solvableDim=sdims[var.label])
    end
  end
  return msg
end

#TODO rename to just prepCliqueMsgDown
"""
    $SIGNATURES

Calculate new and then set the down messages for a clique in Bayes (Junction) tree.
"""
function prepSetCliqueMsgDownConsolidated!( subfg::AbstractDFG,
                                            cliq::TreeClique,
                                            prntDwnMsgs::LikelihoodMessage,
                                            logger=ConsoleLogger();
                                            status::CliqStatus=getCliqueStatus(cliq)  )
  #
  allvars = getCliqVarIdsAll(cliq)
  allprntkeys = collect(keys(prntDwnMsgs.belief))
  passkeys = intersect(allvars, setdiff(allprntkeys,ls(subfg)))
  remainkeys = setdiff(allvars, passkeys)
  newDwnMsgs = LikelihoodMessage(status=status)

  # some msgs are just pass through from parent
  for pk in passkeys
    newDwnMsgs.belief[pk] = prntDwnMsgs.belief[pk]
  end

  # set solvable dimensions
  # sdims = getCliqVariableMoreInitDims(subfg, cliq)

  # other messages must be extracted from subfg
  for mk in remainkeys
    setVari = getVariable(subfg, mk)
    if isInitialized(setVari)
      newDwnMsgs.belief[mk] = TreeBelief(setVari)  #, solvableDim=sdims[mk] )
    end
  end

  # set the downward keys
  with_logger(logger) do
    @info "cliq $(cliq.index), getSetDownMessagesComplete!, allkeys=$(allvars), passkeys=$(passkeys), msgkeys=$(collect(keys(newDwnMsgs.belief)))"
  end

  return newDwnMsgs
end



## =============================================================================
## Multimessage assemblies from multiple cliques 
## =============================================================================
#TODO check as it looks outdated nothing covered after this


#TODO function is currently broken as getUpMsgs does not exist, possibly deprecated?
"""
    $SIGNATURES

Return dictionary of all up belief messages currently in a Bayes `tree`.

Notes
- Returns `::Dict{Int,LikelihoodMessage}`

Related

getUpMsgs
"""
function getTreeCliqUpMsgsAll(tree::AbstractBayesTree)
  allUpMsgs = Dict{Int,LikelihoodMessage}()
  for (idx,cliq) in getCliques(tree)
    msgs = getUpMsgs(cliq)
    allUpMsgs[cliq.index] = LikelihoodMessage()
    for (lbl,msg) in msgs
      # TODO capture the inferred dimension as part of the upward propagation
      allUpMsgs[cliq.index].belief[lbl] = msg
    end
  end
  return allUpMsgs
end

"""
    $SIGNATURES

Convert tree up messages dictionary to a new dictionary relative to variables specific messages and their depth in the tree

Notes
- Used in RoMEPlotting
- Return data in `TempUpMsgPlotting` format:
    Dict{Symbol,   -- is for variable label
      Vector{       -- multiple msgs for the same variable
        Symbol,      -- Clique index
        Int,         -- Depth in tree
        BTD          -- Belief estimate
        inferredDim  -- Information count
      }
"""
function stackCliqUpMsgsByVariable( tree::AbstractBayesTree,
                                    tmpmsgs::Dict{Int, LikelihoodMessage}  )
  #
  # start of the return data structure
  stack = TempUpMsgPlotting()

  # look at all the clique level data
  for (cidx,tmpmsg) in tmpmsgs
    # look at all variables up msg from each clique
    for (sym,msgdim) in tmpmsg.belief
      # create a new object for a particular variable if hasnt been seen before
      if !haskey(stack,sym)
        # FIXME this is an old message type
        stack[sym] = Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}()
      end
      # assemble metadata
      cliq = getCliques(tree,cidx)
      frt = getCliqFrontalVarIds(cliq)[1]
      # add this belief msg and meta data to vector of variable entry
      push!(stack[sym], (frt, getCliqDepth(tree, cliq),msgdim[1], msgdim[2]))
    end
  end

  return stack
end



"""
    $SIGNATURES

Return dictionary of down messages consisting of all frontal and separator beliefs of this clique.

Notes:
- Fetches numerical results from `subdfg` as dictated in `cliq`.
- return LikelihoodMessage
"""
function getCliqDownMsgsAfterDownSolve(subdfg::AbstractDFG, cliq::TreeClique)
  # Dict{Symbol, BallTreeDensity}
  # where the return msgs are contained
  container = LikelihoodMessage() 

  # go through all msgs one by one
  for sym in getCliqAllVarIds(cliq)
    container.belief[sym] = TreeBelief( getVariable(subdfg, sym) )
  end

  # return the result
  return container
end



"""
    $SIGNATURES

Calculate a new down message from the parent.

DevNotes
- FIXME should be handled in CSM
"""
function convertLikelihoodToVector( prntmsgs::Dict{Int, LikelihoodMessage};
                                    logger=SimpleLogger(stdout) )
  #
  # check if any msgs should be multiplied together for the same variable
  
  # FIXME type instability
  msgspervar = Dict{Symbol, Vector{TreeBelief}}()
  for (msgcliqid, msgs) in prntmsgs
    # with_logger(logger) do  #   @info "convertLikelihoodToVector -- msgcliqid=$msgcliqid, msgs.belief=$(collect(keys(msgs.belief)))"  # end
    for (msgsym, msg) in msgs.belief
      # re-initialize with new type
      varType = typeof(msg.softtype)
      # msgspervar = msgspervar !== nothing ? msgspervar : Dict{Symbol, Vector{TreeBelief{varType}}}()
      if !haskey(msgspervar, msgsym)
        # there will be an entire list...
        msgspervar[msgsym] = TreeBelief{varType}[]
      end
      # with_logger(logger) do  @info "convertLikelihoodToVector -- msgcliqid=$(msgcliqid), msgsym $(msgsym), inferdim=$(msg.inferdim)"  # end
      push!(msgspervar[msgsym], msg)
    end
  end

  return msgspervar
end






#
