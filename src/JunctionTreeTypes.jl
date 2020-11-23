

## Bayes Trees

abstract type AbstractBayesTree end

# TODO - see #540 related to indexing and ids
# BayesTree declarations
const BTGdict = GenericIncidenceList{TreeClique,Edge{TreeClique},Array{TreeClique,1},Array{Array{Edge{TreeClique},1},1}}

"""
$(TYPEDEF)

Data structure for the Bayes (Junction) tree, which is used for inference and constructed from a given `::AbstractDFG`.
"""
mutable struct BayesTree <: AbstractBayesTree
  bt::BTGdict
  btid::Int
  cliques::Dict{Int,TreeClique}
  frontals::Dict{Symbol,Int}
  messageChannels::Dict{Int, NamedTuple{(:upMsg, :downMsg),Tuple{Channel{LikelihoodMessage},Channel{LikelihoodMessage}}}}
  variableOrder::Vector{Symbol}
  buildTime::Float64
end

BayesTree() = BayesTree(Graphs.inclist(TreeClique,is_directed=true),
                        0,
                        Dict{Int,TreeClique}(),
                        Dict{AbstractString, Int}(),
                        Dict{Int, NamedTuple{(:upMsg, :downMsg),Tuple{Channel{LikelihoodMessage},Channel{LikelihoodMessage}}}}(),
                        Symbol[],
                        0.0  )
#

#TEMP switch the tree to use NOTE under development don't use MetaBayesTree yet
global UseMetaBayesTree = false
setUseMetaBayesTree(b::Bool) = global UseMetaBayesTree = b
function emptyBayesTree()
  global UseMetaBayesTree
  if UseMetaBayesTree
    @warn "Experimental, do not use yet, MetaBayesTree is under development"
    return MetaBayesTree()
  else
    return BayesTree()
  end
end

# TODO DEV MetaGraphs bayes tree, will potentially also make a LightBayesTree, CloudBayesTree,
"""
$(TYPEDEF)
Data structure for the Bayes (Junction) tree, which is used for inference and constructed from a given `::AbstractDFG`.
"""
mutable struct MetaBayesTree <: AbstractBayesTree
  bt::MetaDiGraph{Int,Float64}
  btid::Int
  # cliques::Dict{Int,TreeClique}
  frontals::Dict{Symbol,Int}
  variableOrder::Vector{Symbol}
  buildTime::Float64
end

MetaBayesTree() = MetaBayesTree(MetaDiGraph{Int,Float64}(), 0, Dict{AbstractString, Int}(), Symbol[], 0.0)

Base.propertynames(x::MetaBayesTree, private::Bool=false) = (:bt, :btid, :cliques, :frontals, :variableOrder, :buildTime)

Base.getproperty(x::MetaBayesTree,f::Symbol) = begin
    if f == :cliques
      if !(@isdefined getCliquesWarnOnce)
        @warn "Maybe don't use cliques field directly, TODO implement add/update/get/delete eg. getClique(tree, cliqId)"
        global getCliquesWarnOnce = true
      end
      d = Dict{Int,Any}()
      for (k,v) in x.bt.vprops
        d[k] = v[:clique]
      end
      return d
    else
      getfield(x,f)
    end
  end

function Base.setproperty!(x::MetaBayesTree, f::Symbol, val)
  if f == :cliques
    if !(@isdefined setCliquesWarnOnce)
      @warn "Maybe don't use cliques field directly, TODO implement add/update/get/delete eg. getClique(tree, cliqId)"
      global setCliquesWarnOnce = true
    end
    for (k,v) in val
      set_prop!(x.bt, k, :clique, v)
    end
  else
    setfield!(x,f,val)
  end
end

function MetaBayesTree(tree::BayesTree)
  # create graph from Graphs.jl adjacency_matrix
  mtree = MetaBayesTree(MetaDiGraph{Int, Float64}(MetaGraphs.SimpleDiGraph(Graphs.adjacency_matrix(tree.bt))), tree.btid, tree.frontals, tree.variableOrder, tree.buildTime)

  #deep copy over properties
  for v in tree.bt.vertices
    # set_prop!(mtree.bt, v.index, :label, deepcopy(v.label))
    set_prop!(mtree.bt, v.index, :clique, deepcopy(v))
  end

  ##  FIXME: Use common location for channels #675 (DF, asking for BTND)
  ##  TODO: placeholder for edge stored Channels
  ## set message passing properties,
  # for e in MetaGraphs.edges(mtree.bt)
  #   set_prop!(mtree.bt, e, :upMsg, Channel{BelieveMessage}(0))
  #   set_prop!(mtree.bt, e, :downMsg, Channel{BelieveMessage}(0))
  # end

  return mtree

end


"""
    $TYPEDEF

Container for upward tree solve / initialization.

DevNotes
- TODO more direct clique access (cliq, parent, children), for multi-process solves
"""
mutable struct CliqStateMachineContainer{BTND, G <: AbstractDFG, InMemG <: InMemoryDFGTypes, BT <: AbstractBayesTree}
  dfg::G
  cliqSubFg::InMemG
  tree::BT
  cliq::TreeClique
  parentCliq::Vector{TreeClique} #TODO deprecate
  childCliqs::Vector{TreeClique} #TODO deprecate
  incremental::Bool
  drawtree::Bool
  dodownsolve::Bool
  delay::Bool
  opts::SolverParams
  refactoring::Dict{Symbol, String}
  oldcliqdata::BTND
  logger::SimpleLogger
  cliqKey::Int
  algorithm::Symbol
end

#TODO use @NamedTuple if julia compat > 1.5
export CSMHistoryTuple
const CSMHistoryTuple =  NamedTuple{(:timestamp, :id, :f, :csmc), Tuple{DateTime, Int, Function, CliqStateMachineContainer}}
const CSMHistory = Vector{CSMHistoryTuple}

Base.show(io::IO, o::CSMHistoryTuple) = print("$(o[1]), $(o[2]), $(o[3])")


function CliqStateMachineContainer( dfg::G,
                                    cliqSubFg::M,
                                    tree::T,
                                    cliq::TreeClique,
                                    parentCliq::Vector{TreeClique},
                                    childCliqs::Vector{TreeClique},
                                    incremental::Bool,
                                    drawtree::Bool,
                                    dodownsolve::Bool,
                                    delay::Bool,
                                    opts::SolverParams,
                                    refactoring::Dict{Symbol,String}=Dict{Symbol,String}(),
                                    oldcliqdata::BTND=BayesTreeNodeData(),
                                    logger::SimpleLogger=SimpleLogger(Base.stdout);
                                    cliqKey::Int = cliq.id,
                                    algoritm::Symbol = :default) where {BTND, G <: AbstractDFG, M <: InMemoryDFGTypes, T <: AbstractBayesTree}
  #
  CliqStateMachineContainer{BTND, G, M, T}( dfg,
                                            cliqSubFg,
                                            tree,
                                            cliq,
                                            parentCliq,
                                            childCliqs,
                                            incremental,
                                            drawtree,
                                            dodownsolve,
                                            delay,
                                            opts,
                                            refactoring,
                                            oldcliqdata,
                                            logger,
                                            cliqKey )
  #
end


function compare(cs1::CliqStateMachineContainer{BTND1, T1, InMemG1, BT1},
                 cs2::CliqStateMachineContainer{BTND2, T2, InMemG2, BT2};
                 skip::Vector{Symbol}=Symbol[] ) where {BTND1, T1 <: AbstractDFG, InMemG1 <: InMemoryDFGTypes, BT1 <: AbstractBayesTree, BTND2, T2 <: AbstractDFG, InMemG2 <: InMemoryDFGTypes, BT2 <: AbstractBayesTree}
  #
  BTND1 == BTND2 ? nothing : @warn("oldcliqdata::$BTND1 != oldcliqdata::$BTND2")
  T1 == T2 ? nothing : @warn("dfg::$T1 != dfg::$T2")
  InMemG1 == InMemG2 ? nothing : @warn("cliqSubFg::$InMemG1 != cliqSubFg::$InMemG2")
  BT1 == BT2 ? nothing : @warn("tree::$BQ1 != tree::$BT2")

  TP = true
  @warn "Skipping compare of CSMC.dfg and .cliqSubFg"
  # TP = TP && compare(cs1.dfg,  cs2.dfg)
  # TP = TP && compare(cs1.cliqSubFg,  cs2.cliqSubFg)
  @warn "Skipping compare of CSMC.tree"
  # TP = TP && compare(cs1.tree,  cs2.tree)
  TP = TP && compare(cs1.cliq,  cs2.cliq)
  TP = TP && compare(cs1.cliqKey,  cs2.cliqKey)
  TP = TP && length(cs1.parentCliq) == length(cs2.parentCliq)
  for i in 1:length(cs1.parentCliq)
    TP = TP && compare(cs1.parentCliq[i],  cs2.parentCliq[i])
  end
  TP = TP && length(cs1.childCliqs) == length(cs2.childCliqs)
  for i in 1:length(cs1.childCliqs)
    TP = TP && compare(cs1.childCliqs[i],  cs2.childCliqs[i])
  end
  TP = TP && compare(cs1.incremental,  cs2.incremental)
  TP = TP && compare(cs1.drawtree,  cs2.drawtree)
  TP = TP && compare(cs1.dodownsolve,  cs2.dodownsolve)
  TP = TP && compare(cs1.delay,  cs2.delay)
  @warn "skipping compare on csmc.opts::SolverParams"
  # TP = TP && compare(cs1.opts,  cs2.opts)
  TP = TP && compare(cs1.refactoring,  cs2.refactoring)
  # TP = TP && compare(cs1.oldcliqdata,  cs2.oldcliqdata)
  # TP = TP && compare(cs1.logger,  cs2.logger)

  return TP
end
