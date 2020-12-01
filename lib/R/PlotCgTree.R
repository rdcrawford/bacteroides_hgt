# ------------------------------------------------------------------------------
# PlotCgTree
# 2020/02/20
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

tree = geneTrees[[8]]
PlotGeneTree( geneTrees[[which( coreGeneClIdxs == boHgtEventData[[1]][[HGT_GENES]][30] )]])
PlotGeneTree = function( tree )
{
  tree = drop.tip( tree, tree$tip.label[ !tree$tip.label %in% genomeNames ] )

  labels = sapply( tree$tip.label, function(x)
  {
    if (  isBx[ which( genomeNames == x ) ] ) return( "red4" )
    return( "darkblue" )
  })
  plot( tree )
  tiplabels( col = labels, pch = 15 ) 
}
for ( i in sample.int( 1384, 10 ) ) PlotGeneTree( geneTrees[[i]] )
PlotCgTree = function( 
  tree, offsetVal, midRoot, plotLegend, algnTips, showLabels, adjVal
  )
{
  # ---- Parse the arguments ---------------------------------------------------
  
  if ( missing( offsetVal ) )  offsetVal  = 0.001
  if ( missing( midRoot ) )    midRoot    = TRUE
  if ( missing( plotLegend ) ) plotLegend = TRUE
  if ( midRoot )               tree       = midpoint.root( tree )
  if ( missing( algnTips ) )
  {
    algnTips = FALSE
  } else {
    algnTips = TRUE
  }
  if ( algnTips ) useLens  = FALSE
  
  if ( missing(showLabels) ) showLabels = TRUE
  
  if ( showLabels ) 
  {
    tipColor = "black"
  } else {
    tipColor = "white"
  }
  
  if ( missing(adjVal) ) adjVal = 25
  
  treePlotChr = 15
  
  # ---- Plot the tree ---------------------------------------------------------
  
  plot.phylo(
    tree, 
    type            = "p", 
    use.edge.length = useLens, 
    font            = 1, 
    label.offset    = offsetVal,
    align.tip.label = algnTips, 
    no.margin       = TRUE, 
    show.tip.label  = TRUE,
    tip.color       = tipColor
    )
  
  # ---- Plot the tip labels ---------------------------------------------------
  
  speciesNames = GetSpeciesLabels( tree$tip.label,   )
  speciesCols = sapply( speciesNames, function(x)
  {
    if ( x == "B. xylanisolvens" ) return( "firebrick" )
    if ( x == "B. ovatus" ) return( "darkblue" )
    return( "grey" )
  })
  
  # Plot the tip labels
  tiplabels(
    pch = treePlotChr,
    col = speciesCols,
    adj = adjVal
    )

  if ( plotLegend )
  {
    legend(
      "bottomright",
      legend = unique(speciesNames),
      col    = unique(speciesCols),
      pch    = 15
      )
  }
}

# ------------------------------------------------------------------------------