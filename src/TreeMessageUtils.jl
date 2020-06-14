# init utils for tree based inference


convert(::Type{BallTreeDensity}, src::TreeBelief) = manikde!(src.val, src.bw[:,1], src.softtype)


## =============================================================================
# helper functions to add tree messages to subgraphs
## =============================================================================


"""
    $SIGNATURES

Modify the `subfg::AbstractDFG` to include `msgs` as priors that are used
during clique inference.

Notes
- May be used initialization or inference, in both upward and downward directions.
- `msgs` are identified by variable label `::Symbol`, and my consist of multiple beliefs.
- Message sets from different cliques are identified by clique id `::Int`.

Related

`deleteMsgFactors!`
"""
function addMsgFactors!(subfg::AbstractDFG,
                        msgs::LikelihoodMessage)
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  for (msym, dm) in msgs.belief
    if msym in svars
      # manifold information is contained in the factor graph DFGVariable object
      fc = addFactor!(subfg, [msym],
              MsgPrior(manikde!(dm.val, dm.bw[:,1], getManifolds(dm.softtype)), dm.inferdim), graphinit=false)
      push!(msgfcts, fc)
    end
  end
  return msgfcts
end


function addMsgFactors!(subfg::AbstractDFG,
                        allmsgs::Dict{Int,LikelihoodMessage} )
  #
  allfcts = DFGFactor[]
  for (cliqid, msgs) in allmsgs
    # do each dict in array separately
    newfcts = addMsgFactors!(subfg, msgs)
    union!( allfcts, newfcts )
  end
  return allfcts
end


# return ::Vector{DFGFactor}
function addMsgFactors_Parametric!(subfg::AbstractDFG,
                                   msgs::LikelihoodMessage)
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
# <<<<<<< Updated upstream
  for (msym, belief) = (msgs.belief)
    if msym in svars
      #TODO covaraince
      #TODO Maybe always use MvNormal
      if size(belief.val)[1] == 1
        msgPrior =  MsgPrior(Normal(belief.val[1], sqrt(belief.bw[1])), belief.inferdim)
      else
        mvnorm = createMvNormal(belief.val[:,1], belief.bw)
        mvnorm == nothing &&
          (return DFGFactor[])
        msgPrior =  MsgPrior(mvnorm, belief.inferdim)
      end
      fc = addFactor!(subfg, [msym], msgPrior, graphinit=false)
      push!(msgfcts, fc)
    end
# =======
#   #convert to Normal or MvNormal
#   varOrder = msgs.cobelief.varlbl
#   Σ = msgs.cobelief.Σ
#   μ = msgs.cobelief.μ
#   if length(μ) == 1
#     dist = Normal(μ[1], sqrt(Σ[1]))
#   else
#     dist = MvNormal(μ, Σ)
#   end
#
#   #TODO add type of factor to message, maybe even send constucted function
#   #TODO inferredDim
#   if length(varOrder) == 1
#     @info "Adding belief msg prior with $dist on $varOrder"
#     fc = addFactor!(subfg, varOrder, MsgPrior(dist, 0.0), graphinit=false)
#   elseif length(varOrder) == 2
#     @error "this only works for linear conditional"
#     fc = addFactor!(subfg, varOrder, LinearConditional(dist), graphinit=false)
#   else
#     error("Oops, not supported")
#   end
#   push!(msgfcts, fc)
#
#   # for (msym, belief) = (msgs.belief)
#   #   if msym in svars
#   #     #TODO covaraince
#   #     #TODO Maybe always use MvNormal
#   #     if size(belief.val)[1] == 1
#   #       msgPrior =  MsgPrior(Normal(belief.val[1], belief.bw[1]), belief.inferdim)
#   #     else
#   #       msgPrior =  MsgPrior(MvNormal(belief.val[:,1], belief.bw), belief.inferdim)
#   #     end
#   #     fc = addFactor!(subfg, [msym], msgPrior, graphinit=false)
#   #     push!(msgfcts, fc)
#   #   end
# >>>>>>> Stashed changes
  end
  return msgfcts
end



"""
    $SIGNATURES

Delete from the subgraph`::AbstractDFG` prior belief `msgs` that could/would be used
during clique inference.

Related

`addMsgFactors!`
"""
function deleteMsgFactors!(subfg::AbstractDFG,
                           fcts::Vector{DFGFactor} )
  #
  for fc in fcts
    deleteFactor!(subfg, fc.label)
  end
end



"""
    $SIGNATURES

Get and return upward belief messages as stored in child cliques from `treel::AbstractBayesTree`.

Notes
- Use last parameter to select the return format.
- Pull model #674

DevNotes
- Consolidate two versions getMsgsUpChildren
"""
function getMsgsUpChildren(fg_::AbstractDFG,
                           treel::AbstractBayesTree,
                           cliq::TreeClique,
                           ::Type{TreeBelief} )
  #
  chld = getChildren(treel, cliq)
  retmsgs = Vector{LikelihoodMessage}(undef, length(chld))
  for i in 1:length(chld)
    retmsgs[i] = getMsgsUpThis(chld[i])
  end
  return retmsgs
end


function getMsgsUpChildren(csmc::CliqStateMachineContainer,
                            ::Type{TreeBelief}=TreeBelief )
  #
  # TODO, replace with single channel stored in csmcs or cliques
  getMsgsUpChildren(csmc.cliqSubFg, csmc.tree, csmc.cliq, TreeBelief)
end


"""
    $SIGNATURES

Get the latest down message from the parent node (without calculating anything).

Notes
- Different from down initialization messages that do calculate new values -- see `prepCliqInitMsgsDown!`.
- Basically converts function `getDwnMsgs` from `Dict{Symbol,BallTreeDensity}` to `Dict{Symbol,Vector{BallTreeDensity}}`.
"""
function getMsgDwnParent(treel::AbstractBayesTree, cliq::TreeClique)
  downmsgs = IntermediateMultiSiblingMessages()
  for prnt in getParent(treel, cliq)
    for (key, bel) in getDwnMsgs(prnt)
      if !haskey(downmsgs, key)
        downmsgs[key] = IntermediateSiblingMessages()
      end
      # TODO insert true inferred dim
      push!(downmsgs[key], bel)
    end
  end
  return downmsgs
end


"""
    $SIGNATURES

Update clique status and notify of the change

Notes
- Assumes users will lock the status state before getting status until after decision whether to update status.
- If so, only unlock after status and condition has been updated.

Dev Notes
- Should be made an atomic transaction
"""
function notifyCliqUpInitStatus!(cliq::TreeClique,
                                 status::Symbol;
                                 logger=ConsoleLogger() )
  #
  cd = getCliqueData(cliq)
  with_logger(logger) do
    tt = split(string(now()), 'T')[end]
    @info "$(tt) $(current_task()), cliq=$(cliq.index), notifyCliqUpInitStatus! -- pre-lock, $(cd.initialized)-->$(status)"
  end
  flush(logger.stream)

  ## TODO only notify if not data structure is not locked by other user (can then remove the hack)
  # Wait until lock can be aquired
  lockUpStatus!(cd)

  cd.initialized = status
  if isready(cd.initUpChannel)
    tkst = take!(cd.initUpChannel)
    # @info "dumping stale cliq=$(cliq.index) status message $(tkst), replacing with $(status)"
  end
  put!(cd.initUpChannel, status)
  cond = getSolveCondition(cliq)
  notify(cond)
    # hack to avoid a race condition  -- remove with atomic lock logic upgrade
    sleep(0.1)
    notify(cond) # getSolveCondition(cliq)

  # TODO unlock
  unlockUpStatus!(cd)
  with_logger(logger) do
    tt = split(string(now()), 'T')[end]
    @info "$(tt) $(current_task()), cliq=$(cliq.index), notifyCliqUpInitStatus! -- unlocked, $(cd.initialized)"
  end

  nothing
end

function notifyCliqDownInitStatus!(cliq::TreeClique,
                                   status::Symbol;
                                   logger=ConsoleLogger() )
  #
  cdat = getCliqueData(cliq)
  with_logger(logger) do
    @info "$(now()) $(current_task()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- pre-lock, new $(cdat.initialized)-->$(status)"
  end

  # take lock for atomic transaction
  lockDwnStatus!(cdat, cliq.index, logger=logger)

  cdat.initialized = status

  if isready(cdat.initDownChannel)
    content = take!(cdat.initDownChannel)
    with_logger(logger) do
      @info "dumping stale cliq=$(cliq.index) status message $(content), replacing with $(status)"
    end
  end
  put!(cdat.initDownChannel, status)
  notify(getSolveCondition(cliq))
    # hack to avoid a race condition
    sleep(0.1)
    notify(getSolveCondition(cliq))

  # unlock for others to proceed
  unlockDwnStatus!(cdat)
  with_logger(logger) do
    @info "$(now()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- unlocked, $(getCliqStatus(cliq))"
  end

  # flush(logger.stream)

  nothing
end




## =============================================================================
## Multimessage assemplies from multiple cliques
## =============================================================================


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
- Return data in `TempUpMsgPlotting` format:
    Dict{Symbol,   -- is for variable label
     Vector{       -- multiple msgs for the same variable
      Symbol,      -- Clique index
      Int,         -- Depth in tree
      BTD          -- Belief estimate
      inferredDim  -- Information count
     }
"""
function stackCliqUpMsgsByVariable(tree::AbstractBayesTree,
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
  container = LikelihoodMessage() # Dict{Symbol,BallTreeDensity}()

  # go through all msgs one by one
  for sym in getCliqAllVarIds(cliq)
    container.belief[sym] = TreeBelief( getVariable(subdfg, sym) )
  end

  # return the result
  return container
end



#
