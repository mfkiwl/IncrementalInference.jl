

"""
    $SIGNATURES

Get the frontal variable IDs `::Int` for a given clique in a Bayes (Junction) tree.
"""
getCliqFrontalVarIds(cliq::Graphs.ExVertex)::Vector{Int} = getData(cliq).frontalIDs

"""
    $(TYPEDSIGNATURES)

Get the frontal variable IDs `::Int` for a given clique in a Bayes (Junction) tree.
"""
getFrontals(cliql::Graphs.ExVertex) = getCliqFrontalVarIds(cliql)



# create a new clique
function addClique!(bt::BayesTree, fg::FactorGraph, varID::Int, condIDs::Array{Int}=Int[])
  bt.btid += 1
  clq = Graphs.add_vertex!(bt.bt, ExVertex(bt.btid,string("Clique",bt.btid)))
  bt.cliques[bt.btid] = clq

  clq.attributes["label"] = ""
  # Specific data container
  setData!(clq, emptyBTNodeData())
  # clq.attributes["data"] = emptyBTNodeData()

  appendClique!(bt, bt.btid, fg, varID, condIDs)
  return clq
end

# generate the label for particular clique -- graphviz drawing
function makeCliqueLabel(fgl::FactorGraph, bt::BayesTree, clqID::Int, api::DataLayerAPI=localapi)
  clq = bt.cliques[clqID]
  flbl = ""
  clbl = ""
  for fr in getData(clq).frontalIDs
    flbl = string(flbl, api.getvertex(fgl,fr).attributes["label"], ",") #fgl.v[fr].
  end
  for cond in getData(clq).conditIDs
    clbl = string(clbl, api.getvertex(fgl,cond).attributes["label"], ",") # fgl.v[cond].
  end
  clq.attributes["label"] = string(flbl, ": ", clbl)
end

# add a conditional ID to clique
function appendConditional(bt::BayesTree, clqID::Int, condIDs::Array{Int,1})
  clq = bt.cliques[clqID]
  getData(clq).conditIDs = union(getData(clq).conditIDs, condIDs)
end

# Add a new frontal variable to clique
function appendClique!(bt::BayesTree, clqID::Int, fg::FactorGraph, varID::Int, condIDs::Array{Int,1}=Int[])
  clq = bt.cliques[clqID]
  var = localapi.getvertex(fg, varID)

  # add frontal variable
  push!(getData(clq).frontalIDs, varID)

  # total dictionary of frontals for easy access
  bt.frontals[var.attributes["label"]] = clqID

  # append to cliq conditionals
  appendConditional(bt, clqID, condIDs)
  makeCliqueLabel(fg, bt, clqID)
  nothing
end


# instantiate a new child clique in the tree
function newChildClique!(bt::BayesTree, fg::FactorGraph, CpID::Int, varID::Int, Sepj::Array{Int,1})
  chclq = addClique!(bt, fg, varID, Sepj)
  parent = bt.cliques[CpID]
  # Staying with Graphs.jl for tree in first stage
  edge = Graphs.make_edge(bt.bt, parent, chclq)
  Graphs.add_edge!(bt.bt, edge)

  return chclq
end

# post order tree traversal and build potential functions
function findCliqueFromFrontal(bt::BayesTree, frtlID::Int)
  for cliqPair in bt.cliques
    id = cliqPair[1]
    cliq = cliqPair[2]
    for frtl in getFrontals(cliq)
      if frtl == frtlID
        return cliq
      end
    end
  end
  error("Clique with desired frontal ID not found")
end


# eliminate a variable for new
function newPotential(tree::BayesTree, fg::FactorGraph, var::Int, prevVar::Int, p::Array{Int,1})
    firvert = localapi.getvertex(fg,var)
    if (length(getData(firvert).separator) == 0)
      if (length(tree.cliques) == 0)
        addClique!(tree, fg, var)
      else
        appendClique!(tree, 1, fg, var) # add to root
      end
    else
      Sj = getData(firvert).separator
      # find parent clique Cp that containts the first eliminated variable of Sj as frontal
      firstelim = 99999999999
      for s in Sj
        temp = something(findfirst(isequal(s), p), 0) # findfirst(p, s)
        if (temp < firstelim)
          firstelim = temp
        end
      end
      felbl = localapi.getvertex(fg, p[firstelim]).attributes["label"]
      CpID = tree.frontals[felbl]
      # look to add this conditional to the tree
      unFC = union(getCliqFrontalVarIds(tree.cliques[CpID]), getCliqSeparatorVarIds(tree.cliques[CpID]))
      if (sort(unFC) == sort(Sj))
        appendClique!(tree, CpID, fg, var)
      else
        newChildClique!(tree, fg, CpID, var, Sj)
      end
    end
end

# build the whole tree in batch format
function buildTree!(tree::BayesTree, fg::FactorGraph, p::Array{Int,1})
  rp = reverse(p,dims=1) # flipdim(p, 1)
  prevVar = 0
  for var in rp
    newPotential(tree, fg, var, prevVar, p)
    prevVar = var
  end
end

function showTree(;filepath::String="/tmp/bt.pdf",
                   viewerapp::String="evince"  )
  #
  try
    @async run(`$(viewerapp) $(filepath)`)
  catch ex
    @warn "not able to show via $(viewerapp) $(filepath)"
    @show ex
    @show stacktrace()
  end
end

function drawTree(treel::BayesTree;
                  show::Bool=false,                  # must remain false for stability and automated use in solver
                  filepath::String="/tmp/bt.pdf",
                  viewerapp::String="evince",
                  imgs::Bool=false )
  #
  fext = split(filepath, '.')[end]
  fpwoext = split(filepath, '.')[end-1]

  # modify a deepcopy
  btc = deepcopy(treel)
  for (cid, cliq) in btc.cliques
    if imgs
      firstlabel = split(cliq.attributes["label"],',')[1]
      spyCliqMat(cliq, suppressprint=true) |> exportimg("/tmp/$firstlabel.png")
      cliq.attributes["image"] = "/tmp/$firstlabel.png"
      cliq.attributes["label"] = ""
    end
    delete!(cliq.attributes, "data")
  end

  fid = IOStream("")
  try
    fid = open("$(fpwoext).dot","w+")
    write(fid,to_dot(btc.bt))
    close(fid)
    run(`dot $(fpwoext).dot -T$(fext) -o $(filepath)`)
  catch ex
    @warn ex
    @show stacktrace()
  finally
    close(fid)
  end

  show ? showTree(viewerapp=viewerapp, filepath=filepath) : nothing
end



## Find batch belief propagation solution
function prepBatchTree!(fg::FactorGraph;
                        ordering::Symbol=:qr,
                        drawpdf::Bool=false,
                        show::Bool=false,
                        filepath::String="/tmp/bt.pdf",
                        viewerapp::String="evince",
                        imgs::Bool=false  )
  #
  p = IncrementalInference.getEliminationOrder(fg, ordering=ordering)
  println()
  fge = deepcopy(fg)
  println("Building Bayes net...")
  buildBayesNet!(fge, p)

  tree = emptyBayesTree()
  buildTree!(tree, fge, p)

  # println("Bayes Net")
  # sleep(0.1)
  #fid = open("bn.dot","w+")
  #write(fid,to_dot(fge.bn))
  #close(fid)

  println("Find potential functions for each clique")
  cliq = tree.cliques[1] # start at the root
  buildCliquePotentials(fg, tree, cliq); # fg does not have the marginals as fge does

  # now update all factor graph vertices used for this tree
  for (id,v) in fg.g.vertices
    dlapi.updatevertex!(fg, v)
  end

  # GraphViz.Graph(to_dot(tree.bt))
  # Michael reference -- x2->x1, x2->x3, x2->x4, x2->l1, x4->x3, l1->x3, l1->x4
  #Michael reference 3sig -- x2l1x4x3    x1|x2
  println("Bayes Tree")
  if drawpdf
    drawTree(tree, show=show, filepath=filepath, viewerapp=viewerapp, imgs=imgs)
  end

  return tree
end

function resetData!(vdata::VariableNodeData)::Nothing
  vdata.eliminated = false
  vdata.BayesNetOutVertIDs = Int[]
  vdata.BayesNetVertID = 0
  vdata.separator = Int[]
  nothing
end

function resetData!(vdata::FunctionNodeData)::Nothing
  vdata.eliminated = false
  vdata.potentialused = false
  nothing
end

function resetFactorGraphNewTree!(fgl::FactorGraph)::Nothing
  for (id, v) in fgl.g.vertices
    resetData!(getData(v))
    localapi.updatevertex!(fgl, v)
  end
  nothing
end

"""
    $(SIGNATURES)

Build a completely new Bayes (Junction) tree, after first wiping clean all
temporary state in fg from a possibly pre-existing tree.
"""
function wipeBuildNewTree!(fg::FactorGraph;
                           ordering::Symbol=:qr,
                           drawpdf::Bool=false,
                           show::Bool=false,
                           filepath::String="/tmp/bt.pdf",
                           viewerapp::String="evince",
                           imgs::Bool=false  )::BayesTree
  #
  resetFactorGraphNewTree!(fg);
  return prepBatchTree!(fg, ordering=ordering, drawpdf=drawpdf, show=show, filepath=filepath, viewerapp=viewerapp, imgs=imgs);
end

"""
    $(SIGNATURES)

Return the Graphs.ExVertex node object that represents a clique in the Bayes
(Junction) tree, as defined by one of the frontal variables `frt<:AbstractString`.
"""
function whichCliq(bt::BayesTree, frt::T) where {T <: AbstractString}
  bt.cliques[bt.frontals[frt]]
end
whichCliq(bt::BayesTree, frt::Symbol) = whichCliq(bt, string(frt))

"""
    $SIGNATURES

Return the Graphs.ExVertex node object that represents a clique in the Bayes
(Junction) tree, as defined by one of the frontal variables `frt::Symbol`.
"""
getCliq(bt::BayesTree, frt::Symbol) = whichCliq(bt, string(frt))



"""
    $(SIGNATURES)

Set the upward passing message for Bayes (Junction) tree clique `cliql`.
"""
function setUpMsg!(cliql::ExVertex, msgs::Dict{Symbol, BallTreeDensity})
  getData(cliql).upMsg = msgs
end

"""
    $(SIGNATURES)

Set the downward passing message for Bayes (Junction) tree clique `cliql`.
"""
function setDwnMsg!(cliql::ExVertex, msgs::Dict{Symbol, BallTreeDensity})
  getData(cliql).dwnMsg = msgs
end

"""
    $(SIGNATURES)

Return the last up message stored in `cliq` of Bayes (Junction) tree.
"""
function upMsg(cliq::Graphs.ExVertex)
  getData(cliq).upMsg
end
function upMsg(btl::BayesTree, sym::Symbol)
  upMsg(whichCliq(btl, sym))
end

"""
    $(SIGNATURES)

Return the last up message stored in `cliq` of Bayes (Junction) tree.
"""
getUpMsgs(btl::BayesTree, sym::Symbol) = upMsg(btl, sym)
getUpMsgs(cliql::Graphs.ExVertex) = upMsg(cliql)



"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
function dwnMsg(cliq::Graphs.ExVertex)
  getData(cliq).dwnMsg
end
function dwnMsg(btl::BayesTree, sym::Symbol)
  upMsg(whichCliq(btl, sym))
end

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
getDwnMsgs(btl::BayesTree, sym::Symbol) = dwnMsg(btl, sym)
getDwnMsgs(cliql::Graphs.ExVertex) = dwnMsg(cliql)


function appendUseFcts!(usefcts, lblid::Int, fct::Graphs.ExVertex, fid::Int)
  for tp in usefcts
    if tp == fct.index
      return nothing
    end
  end
  tpl = fct.index
  push!(usefcts, tpl )
  nothing
end

"""
    $SIGNATURES

Return list of factors which depend only on variables in variable list in factor
graph -- i.e. among variables.

Notes
-----
* `unused::Bool=true` will disregard factors already used -- i.e. disregard where `potentialused=true`
"""
function getFactorsAmongVariablesOnly(fgl::FactorGraph,
                                      varlist::Vector{Symbol};
                                      unused::Bool=true  )
  # collect all factors attached to variables
  prefcts = Symbol[]
  for var in varlist
    union!(prefcts, ls(fgl, var))
  end

  almostfcts = Symbol[]
  if unused
    # now check if those factors have already been added
    for fct in prefcts
      vert = getVert(fgl, fct, nt=:fct)
      if !getData(vert).potentialused
        push!(almostfcts, fct)
      end
    end
  else
    almostfcts = prefcts
  end

  # Select factors that have all variables in this clique var list
  usefcts = Symbol[]
  for fct in almostfcts
    if length(setdiff(lsf(fgl, fct), varlist)) == 0
      push!(usefcts, fct)
    end
  end

  return usefcts
end

# function getCliqFactorsFromFrontals(fgl::FactorGraph,
#                                     varlist::Vector{Symbol};
#                                     unused::Bool=true,
#                                     api::DataLayerAPI=localapi )
#   #
#   frtl = getData(cliq).frontalIDs
#   cond = getData(cliq).conditIDs
#   allids = [frtl;cond]
#
#   usefcts = Int[]
#   for fid in frtl
#     usefcts = Int[]
#     for fct in api.outneighbors(fgl, getVert(fg,fid, api=api))
# #         if getData(fct).potentialused!=tfrue
# #             loutn = localapi.outneighbors(fg, fct)
# #             if length(loutn)==1
# #                 appendUseFcts!(usefcts, fg.IDs[Symbol(loutn[1].label)], fct, fid)
# #                 fct.attributes["data"].potentialused = true
# #                 localapi.updatevertex!(fg, fct)
# #             end
# #             for sepSearch in loutn
# #                 sslbl = Symbol(sepSearch.label)
# #                 if (fg.IDs[sslbl] == fid)
# #                     continue # skip the fid itself
# #                 end
# #                 sea = findmin(abs.(allids .- fg.IDs[sslbl]))
# #                 if sea[1]==0.0
# #                     appendUseFcts!(usefcts, fg.IDs[sslbl], fct, fid)
# #                     fct.attributes["data"].potentialused = true
# #                     localapi.updatevertex!(fg, fct)
# #                 end
# #             end
# #         end
#     end
# #     getData(cliq).potentials = union(getData(cliq).potentials, usefcts)
#   end
#
#   return usefcts
# end

"""
    $SIGNATURES

Return `::Bool` on whether factor `fct<:FunctorInferenceType` is a partial constraint.
"""
isPartial(fcf::T) where {T <: FunctorInferenceType} = :partial in fieldnames(T)
function isPartial(fct::Graphs.ExVertex)
  fcf = getFactor(fct)
  isPartial(fcf)
end

"""
    $SIGNATURES

Get and set the potentials for a particular `cliq` in the Bayes (Junction) tree.
"""
function getCliquePotentials!(fg::FactorGraph,
                              bt::BayesTree,
                              cliq::Graphs.ExVertex  )
  #
  frtl = getData(cliq).frontalIDs
  cond = getData(cliq).conditIDs
  allids = [frtl;cond]

  if true
    varlist = Symbol[]
    for id in allids
      push!(varlist, getSym(fg, id))
    end
    fctsyms = getFactorsAmongVariablesOnly(fg, varlist, unused=true )
    for fsym in fctsyms
      push!(getData(cliq).potentials, fg.fIDs[fsym])
      fct = getVert(fg, fsym, nt=:fct)
      getData(fct).potentialused = true
      push!(getData(cliq).partialpotential, isPartial(fct))
    end
  # else
  #   frontal only factor approach
  # end

  end

  # check if any of the factors are partial constraints
  # prfctsids = getCliqFactorIds(cliq)[prfcts]
  # len = length(prfctsids)
  # fulls = zeros(Bool, len)
  # for i in 1:len
  #   fulls = !isPartial()
  # end


  nothing
end

function getCliquePotentials!(fg::FactorGraph, bt::BayesTree, chkcliq::Int)
    getCliquePotentials!(fg, bt.cliques[chkcliq])
end

function cliqPotentialIDs(cliq::Graphs.ExVertex)
  potIDs = Int[]
  for idfct in getData(cliq).potentials
    push!(potIDs,idfct)
  end
  return potIDs
end

function collectSeparators(bt::BayesTree, cliq::Graphs.ExVertex)
  allseps = Int[]
  for child in out_neighbors(cliq, bt.bt)#tree
      allseps = [allseps; getData(child).conditIDs]
  end
  return allseps
end

"""
    $SIGNATURES

Return boolean matrix of factor by variable (row by column) associations within
clique, corresponds to order presented by `getCliqFactorIds` and `getCliqAllVarIds`.
"""
function getCliqAssocMat(cliq::Graphs.ExVertex)
  getData(cliq).cliqAssocMat
end

"""
    $SIGNATURES

Return boolean matrix of upward message singletons (i.e. marginal priors) from
child cliques.  Variable order corresponds to `getCliqAllVarIds`.
"""
function getCliqMsgMat(cliq::Graphs.ExVertex)
  getData(cliq).cliqMsgMat
end

"""
    $SIGNATURES

Return boolean matrix of factor variable associations for a clique, optionally
including (`showmsg::Bool=true`) the upward message singletons.  Variable order
corresponds to `getCliqAllVarIds`.
"""
function getCliqMat(cliq::Graphs.ExVertex; showmsg::Bool=true)
  assocMat = getCliqAssocMat(cliq)
  msgMat = getCliqMsgMat(cliq)
  mat = showmsg ? [assocMat;msgMat] : assocMat
  return mat
end


"""
    $SIGNATURES

Get `cliq` separator (a.k.a. conditional) variable ids`::Int`.
"""
function getCliqSeparatorVarIds(cliq::Graphs.ExVertex)::Vector{Int}
  getData(cliq).conditIDs
end

"""
    $SIGNATURES

Get `cliq` potentials (factors) ids`::Int`.
"""
function getCliqFactorIds(cliq::Graphs.ExVertex)::Vector{Int}
  getData(cliq).potentials
end

"""
    $SIGNATURES

Get all `cliq` variable ids`::Int`.
"""
function getCliqAllVarIds(cliq::Graphs.ExVertex)::Vector{Int}
  frtl = getCliqFrontalVarIds(cliq)
  cond = getCliqSeparatorVarIds(cliq)
  [frtl;cond]
end

"""
    $SIGNATURES

Get variable ids`::Int` with prior factors associated with this `cliq`.
"""
function getCliqVarIdsPriors(cliq::Graphs.ExVertex,
                             allids::Vector{Int}=getCliqAllVarIds(cliq),
                             partials::Bool=true  )::Vector{Int}
  # get ids with prior factors associated with this cliq
  amat = getCliqAssocMat(cliq)
  prfcts = sum(amat, dims=2) .== 1

  # remove partials as per request
  !partials ? nothing : (prfcts .&= getData(cliq).partialpotential)

  # return variable ids in `mask`
  mask = sum(amat[prfcts[:],:], dims=1)[:] .> 0
  return allids[mask]
end

"""
    $SIGNATURES

Return the number of factors associated with each variable in `cliq`.
"""
getCliqNumAssocFactorsPerVar(cliq::Graphs.ExVertex)::Vector{Int} = sum(getCliqAssocMat(cliq), dims=1)[:]



"""
    $SIGNATURES

Get `cliq` variable IDs with singleton factors -- i.e. both in clique priors and up messages.
"""
function getCliqVarSingletons(cliq::Graphs.ExVertex,
                              allids::Vector{Int}=getCliqAllVarIds(cliq),
                              partials::Bool=true  )::Vector{Int}
  # get incoming upward messages (known singletons)
  mask = sum(getCliqMsgMat(cliq),dims=1)[:] .>= 1
  upmsgids = allids[mask]

  # get ids with prior factors associated with this cliq
  prids = getCliqVarIdsPriors(cliq, getCliqAllVarIds(cliq), partials)

  # return union of both lists
  return union(upmsgids, prids)
end

"""
    $SIGNATURES

Based on a push model from child cliques that should have already completed their computation.
"""
function getCliqInitUpMsgs(cliq::Graphs.ExVertex)
  getData(cliq).upInitMsg
end

function setCliqUpInitMsg!(cliq::Graphs.ExVertex, childid::Int, msg::Dict{})
  getData(cliq).upInitMsg[childid] = msg
end

function isCliqInitialized(cliq::Graphs.ExVertex)::Bool
  return getData(cliq).initialized in [:initialized; :upsolved]
end

"""
    $SIGNATURES

Return the most likely  ordering for initializing factor (assuming up solve
sequence).
"""
function getCliqInitVarOrderUp(cliq::Graphs.ExVertex)
  # rules to explore dimension from one to the other?

  # get all variable ids and number of associated factors
  allids = getCliqAllVarIds(cliq)
  nfcts = getCliqNumAssocFactorsPerVar(cliq)

  # get priors and singleton message variables (without partials)
  prids = getCliqVarIdsPriors(cliq, getCliqAllVarIds(cliq), false)

  # get current up msgs in the init process (now have all singletons)
  upmsgs = getCliqInitUpMsgs(cliq)
  upmsgids = collect(keys(upmsgs))

  # all singleton variables
  singids = union(prids, upmsgids)

  # add msg marginal prior (singletons) to number of factors
  for msid in upmsgids
    nfcts[msid .== allids] .+= 1
  end

  # sort permutation order for increasing number of factor association
  nfctsp = sortperm(nfcts)
  sortedids = allids[nfctsp]

  # organize the prior variables separately with asceding factor count
  initorder = zeros(Int, 0)
  for id in sortedids
    if id in singids
      push!(initorder, id)
    end
  end
  # in ascending order of number of factors
  for id in sortedids
    if !(id in initorder)
      push!(initorder, id)
    end
  end
  return initorder
end

"""
    $SIGNATURES

Return true if clique has completed the local upward direction inference procedure.
"""
isUpInferenceComplete(cliq::Graphs.ExVertex) = getData(cliq).upsolved

function areCliqVariablesInitialized(fgl::FactorGraph, cliq::Graphs.ExVertex)
  allids = getCliqAllVarIds(cliq)
  isallinit = true
  for vid in allids
    var = getVert(fgl, vid, api=localapi)
    isallinit &= isInitialized(var)
  end
  isallinit
end

"""
   $SIGNATURES

Determine if this `cliq` has been fully initialized and child cliques have completed their full upward inference.
"""
function isCliqReadyInferenceUp(fgl::FactorGraph, tree::BayesTree, cliq::Graphs.ExVertex)
  isallinit = areCliqVariablesInitialized(fgl, cliq)

  # check that all child cliques have also completed full up inference.
  for chl in getChildren(tree, cliq)
    @show isallinit &= isUpInferenceComplete(chl)
  end
  return isallinit
end



"""
    $SIGNATURES

Perform cliq initalization calculation based on current state of the tree and factor graph,
using upward message passing logic.

> NOTE WORK IN PROGRESS

Notes
- Return either of (:initialized, :upsolved, :needdownmsg, :badinit)
"""
function doCliqAutoInitUp!(fgl::FactorGraph,
                           tree::BayesTree,
                           cliq::Graphs.ExVertex;
                           up_solve_if_able::Bool=true  )::Symbol
  #
  # init up msg has special procedure for incomplete messages
  retmsg = :badinit
  msg = Dict{Symbol, BallTreeDensity}()

  # structure for all up message densities computed during this initialization procedure.
  varorder = getCliqInitVarOrderUp(cliq)

  # do physical inits
  count = 1
  while count > 0
    count = 0
    for vid in varorder
      var = getVert(fgl, vid, api=localapi)
      isinit = isInitialized(var)
      doautoinit!(fgl, ExVertex[var;], api=localapi)
      isinit == isInitialized(var) ? nothing : (count += 1)
    end
  end

  # next default return type
  retmsg = :needdownmsg

  # check if all cliq vars have been initialized so that full inference can occur on clique
  isinit = areCliqVariablesInitialized(fgl, cliq)
  # might fail while waiting for other cliques to initialize.
  if isinit
    retmsg = :initialized
    if up_solve_if_able
      csym = Symbol(getVert(fgl,getCliqFrontalVarIds(cliq)[1],api=localapi).label)
      approxCliqMarginalUp!(fgl, tree, csym, false)
      retmsg = :solved
    end
  end

  # construct init msg to place in parent as initialized separator variables
  for vid in getCliqSeparatorVarIds(cliq)
    var = getVert(fgl, vid, api=localapi)
    if isInitialized(var)
      msg[Symbol(var.label)] = getKDE(var)
    end
  end

  # put the init result in the parent cliq.
  prnt = getParent(tree, cliq)
  # not a root clique
  if length(prnt) > 0
    setCliqUpInitMsg!(prnt[1], cliq.index, msg)
  end

  getData(cliq).initialized = retmsg
  return retmsg
end



# """
#     $SIGNATURES
#
# Based on current status in factor graph, determine if initialization of requested
# variable is possible.
# """
# function getCliqInitVarPossible(cliq::Graphs.ExVertex, varid::Int)
#
#   factorCanInitFromOtherVars(cliq, fctid)
#
# end

# """
#     $SIGNATURES
#
# Determine if a clique Chapman-Kolmogorov computation can be achieved,
# alongide additional message singletons that might be available from caller.
# """
# function calcCliqTotalSolvePossible(cliq::Graphs.ExVertex;
#                                     allids::Vector{Int}=getCliqAllVarIds(cliq),
#                                     availablemsgs::Vector{Bool}=zeros(Bool,length(allids)) )::Tuple{Bool, Vector{Bool}}
#   # return list of all initable variables in cliq (default is false)
#   initable = zeros(Bool, length(allids))
#
#   # what is the initialization order
#   initorder = getCliqInitVarOrderUp(cliq::Graphs.ExVertex)
#
#   # check if all variables can be initialized
#   for i in 1:length(allids)
#     if allids[i] in initorder
#       if getCliqInitVarPossible(cliq, allids[i])
#         initable[i] = true
#       end
#     end
#   end
#
#   # would be fully initializable if all initable are true
#   return sum(initable)==length(allids), initable
# end


function compCliqAssocMatrices!(fgl::FactorGraph, bt::BayesTree, cliq::Graphs.ExVertex)
  frtl = getCliqFrontalVarIds(cliq)
  cond = getCliqSeparatorVarIds(cliq)
  inmsgIDs = collectSeparators(bt, cliq)
  potIDs = cliqPotentialIDs(cliq)
  # Construct associations matrix here
  # matrix has variables are columns, and messages/constraints as rows
  cols = [frtl;cond]
  getData(cliq).inmsgIDs = inmsgIDs
  getData(cliq).potIDs = potIDs
  cliqAssocMat = Array{Bool,2}(undef, length(potIDs), length(cols))
  cliqMsgMat = Array{Bool,2}(undef, length(inmsgIDs), length(cols))
  fill!(cliqAssocMat, false)
  fill!(cliqMsgMat, false)
  for j in 1:length(cols)
    for i in 1:length(inmsgIDs)
      if cols[j] == inmsgIDs[i]
        cliqMsgMat[i,j] = true
      end
    end
    for i in 1:length(potIDs)
      idfct = getData(cliq).potentials[i]
      if idfct == potIDs[i] # sanity check on clique potentials ordering
        for vertidx in getData(getVert(fgl, idfct, api=localapi)).fncargvID
        # for vertidx in getData(getVertNode(fgl, idfct)).fncargvID
          if vertidx == cols[j]
            cliqAssocMat[i,j] = true
          end
        end
      else
        prtslperr("compCliqAssocMatrices! -- potential ID ordering was lost")
      end
    end
  end
  getData(cliq).cliqAssocMat = cliqAssocMat
  getData(cliq).cliqMsgMat = cliqMsgMat
  nothing
end


function countSkips(bt::BayesTree)
  skps = 0
  for cliq in bt.cliques
    m = getCliqMat(cliq[2])
    mi = map(Int,m)
    skps += sum(map(Int,sum(mi, dims=1) .== 1))
  end
  return skps
end

function skipThroughMsgsIDs(cliq::Graphs.ExVertex)
  cliqdata = getData(cliq)
  numfrtl1 = floor(Int,length(cliqdata.frontalIDs)+1)
  condAssocMat = cliqdata.cliqAssocMat[:,numfrtl1:end]
  condMsgMat = cliqdata.cliqMsgMat[:,numfrtl1:end]
  mat = [condAssocMat;condMsgMat];
  mab = sum(map(Int,mat),dims=1) .== 1
  mabM = sum(map(Int,condMsgMat),dims=1) .== 1
  mab = mab .& mabM
  # rang = 1:size(condMsgMat,2)
  msgidx = cliqdata.conditIDs[vec(collect(mab))]
  return msgidx
end

function directPriorMsgIDs(cliq::Graphs.ExVertex)
  frtl = getData(cliq).frontalIDs
  cond = getData(cliq).conditIDs
  cols = [frtl;cond]
  mat = getCliqMat(cliq, showmsg=true)
  singr = sum(map(Int,mat),dims=2) .== 1
  rerows = collect(1:length(singr))
  b = vec(collect(singr))
  rerows2 = rerows[b]
  sumsrAc = sum(map(Int,mat[rerows2,:]),dims=1)
  sumc = sum(map(Int,mat),dims=1)
  pmSkipCols = (sumsrAc - sumc) .== 0
  return cols[vec(collect(pmSkipCols))]
end

function directFrtlMsgIDs(cliq::Graphs.ExVertex)
  numfrtl = length(getData(cliq).frontalIDs)
  frntAssocMat = getData(cliq).cliqAssocMat[:,1:numfrtl]
  frtlMsgMat = getData(cliq).cliqMsgMat[:,1:numfrtl]
  mat = [frntAssocMat; frtlMsgMat];
  mab = sum(map(Int,mat),dims=1) .== 1
  mabM = sum(map(Int,frtlMsgMat),dims=1) .== 1
  mab = mab .& mabM
  return getData(cliq).frontalIDs[vec(collect(mab))]
end

function directAssignmentIDs(cliq::Graphs.ExVertex)
  # NOTE -- old version been included in iterated variable stack
  assocMat = getData(cliq).cliqAssocMat
  msgMat = getData(cliq).cliqMsgMat
  mat = [assocMat;msgMat];
  mab = sum(map(Int,mat),dims=1) .== 1
  mabA = sum(map(Int,assocMat),dims=1) .== 1
  mab = mab .& mabA
  frtl = getData(cliq).frontalIDs
  cond = getData(cliq).conditIDs
  cols = [frtl;cond]
  return cols[vec(collect(mab))]
  # also calculate how which are conditionals
end

function mcmcIterationIDs(cliq::Graphs.ExVertex)
  mat = getCliqMat(cliq)
  # assocMat = getData(cliq).cliqAssocMat
  # msgMat = getData(cliq).cliqMsgMat
  # mat = [assocMat;msgMat];

  sum(sum(map(Int,mat),dims=1)) == 0 ? error("mcmcIterationIDs -- unaccounted variables") : nothing
  mab = 1 .< sum(map(Int,mat),dims=1)
  cols = getCliqAllVarIds(cliq)

  # must also include "direct variables" connected through projection only
  directvars = directAssignmentIDs(cliq)
  usset = union(directvars, cols[vec(collect(mab))])
  # NOTE -- fix direct vs itervar issue, DirectVarIDs against Iters should also Iter
  # NOTE -- using direct then mcmcIter ordering to prioritize non-msg vars first
  return setdiff(usset, getData(cliq).directPriorMsgIDs)
end

function getCliqMatVarIdx(cliq::Graphs.ExVertex, varid::Int, allids=getCliqAllVarIds(cliq) )
  len = length(allids)
  [1:len;][allids .== varid][1]
end

"""
    $SIGNATURES

Determine and return order list of variable ids required for minibatch Gibbs iteration inside `cliq`.

Notes
* Singleton factors (priors and up messages) back of the list
* least number of associated factor variables earlier in list
"""
function mcmcIterationIdsOrdered(cliq::Graphs.ExVertex)
  # get unordered iter list
  alliter = mcmcIterationIDs(cliq)

  # get all singletons
  allsings = getCliqVarSingletons(cliq)
  singletonvars = intersect(alliter, allsings)

  # get all non-singleton iters
  nonsinglvars = setdiff(alliter, singletonvars)

  # sort nonsingletons ascending number of factors
  mat = getCliqMat(cliq)
  lenfcts = sum(mat, dims=1)
  nonslen = zeros(length(nonsinglvars))
  for i in 1:length(nonsinglvars)
    varid = nonsinglvars[i]
    varidx = getCliqMatVarIdx(cliq, varid)
    nonslen[i] = lenfcts[varidx]
  end
  p = sortperm(nonslen)
  ascnons = nonsinglvars[p]

  # sort singleton vars ascending number of factors
  singslen = zeros(length(singletonvars))
  for i in 1:length(singletonvars)
    varid = singletonvars[i]
    varidx = getCliqMatVarIdx(cliq, varid)
    singslen[i] = lenfcts[varidx]
  end
  p = sortperm(singslen)
  ascsing = singletonvars[p]

  return [ascnons; ascsing]
end

"""
    $(SIGNATURES)

Prepare the variable IDs for nested clique Gibbs mini-batch calculations, by assembing these clique data fields:
- `directPriorMsgIDs`
- `directvarIDs`
- `itervarIDs`
- `msgskipIDs`
- `directFrtlMsgIDs`

"""
function setCliqMCIDs!(cliq::Graphs.ExVertex)
  getData(cliq).directPriorMsgIDs = directPriorMsgIDs(cliq)

  # NOTE -- directvarIDs are combined into itervarIDs
  getData(cliq).directvarIDs = directAssignmentIDs(cliq)
  # TODO find itervarIDs that have upward child singleton messages and update them last in iter list
  getData(cliq).itervarIDs = mcmcIterationIdsOrdered(cliq)  #mcmcIterationIDs(cliq)

  getData(cliq).msgskipIDs = skipThroughMsgsIDs(cliq)
  getData(cliq).directFrtlMsgIDs = directFrtlMsgIDs(cliq)

  # TODO add initialization sequence var id list too

  nothing
end


# post order tree traversal and build potential functions
function buildCliquePotentials(fg::FactorGraph, bt::BayesTree, cliq::Graphs.ExVertex)
    for child in out_neighbors(cliq, bt.bt)#tree
        buildCliquePotentials(fg, bt, child)
    end
    @info "Get potentials $(cliq.attributes["label"])"
    getCliquePotentials!(fg, bt, cliq);

    compCliqAssocMatrices!(fg, bt, cliq);
    setCliqMCIDs!(cliq);

    nothing
end

"""
    $(SIGNATURES)

Return a vector of child cliques to `cliq`.
"""
function childCliqs(treel::BayesTree, cliq::Graphs.ExVertex)
    childcliqs = Vector{Graphs.ExVertex}(undef, 0)
    for cl in Graphs.out_neighbors(cliq, treel.bt)
        push!(childcliqs, cl)
    end
    return childcliqs
end
function childCliqs(treel::BayesTree, frtsym::Symbol)
  childCliqs(treel,  whichCliq(treel, frtsym))
end

"""
    $(SIGNATURES)

Return a vector of child cliques to `cliq`.
"""
getChildren(treel::BayesTree, frtsym::Symbol) = childCliqs(treel, frtsym)
getChildren(treel::BayesTree, cliq::Graphs.ExVertex) = childCliqs(treel, cliq)

"""
    $(SIGNATURES)

Return `cliq`'s parent clique.
"""
function parentCliq(treel::BayesTree, cliq::Graphs.ExVertex)
    Graphs.in_neighbors(cliq, treel.bt)
end
function parentCliq(treel::BayesTree, frtsym::Symbol)
  parentCliq(treel,  whichCliq(treel, frtsym))
end

"""
    $(SIGNATURES)

Return `cliq`'s parent clique.
"""
getParent(treel::BayesTree, afrontal::Union{Symbol, Graphs.ExVertex}) = parentCliq(treel, afrontal)
