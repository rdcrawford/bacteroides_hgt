# ------------------------------------------------------------------------------
# 
# 2020/07/27
# Ryan D. Crawford
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( Rcpp )
library( cognac )
library( ape )
library( maps )
library( phytools )
library( RColorBrewer )

# ---- Constant Declarations ---------------------------------------------------

HGT_THRESH    = 0.1
G_ANNOT       = 2
G_START       = 4
G_END         = 5
G_CLUST       = 7
MAX_GENE_DIST = 5000
GENOME_NAME   = 1
CLUST_IDX     = 2
HGT_GENES     = 3
CG_DISTS      = 4
ACC_GENES     = 5
ACC_GENE_PRS  = 6
HGT_ACC_GENES = 7
GENE_ANNOTS   = 8

# ------------------------------------------------------------------------------

# Load the data on the gene aligments made with cognac 
load( "data/2020_09_23_cognac_run_data.Rdata" )

rData = "2020_07_30_analyze_tree_order.Rdata"

# ------------------------------------------------------------------------------

bxCols = 
  colorRampPalette( c( "indianred", "red4" ) )( length(speciesClusts[[1]]))
boCols = colorRampPalette( c( "cyan", "navy" ) )( length(speciesClusts[[2]]))
names( bxCols ) = speciesClusts[[1]]
names( boCols ) = speciesClusts[[2]]

genomeIds = c( speciesClusts[[1]], speciesClusts[[2]] )

sortedMat = t( sapply( 1:length(geneTrees),
  function(i) GetTipLabelOrder( geneTrees[[i]], genomeIds )
  ) )

pheatmap::pheatmap(
  sortedMat,
  color = c( bxCols, boCols )
  )

gff = genePtr$gfList[[ 16 ]]
geneOrder = sapply( 1:length( genePtr$clustList ),
  function(i) which( gff[ , 1 ] %in% genePtr$clustList[[i]] )
  )
sortedMat = sortedMat[ order( geneOrder ), ] 

sortedGenes = geneOrder[ order( geneOrder ) ]

# ------------------------------------------------------------------------------

sGenes = geneOrder[ order( geneOrder ) ]

pheatmap::pheatmap(
  sortedMat,
  color = c( bxCols, boCols ),
  cluster_rows = FALSE
  )

gapPos = 1
i = 1
while ( i < length( sortedGenes ) - 1)
{
  if ( sortedGenes[ i + 1 ] - sortedGenes[ i ] != 1 ) 
    gapPos[ length(gapPos) + 1] = i
  i = i + 1
}

numGenes = sapply( 1:( length(gapPos) - 1 ), 
  function(i) gapPos[ i + 1 ] - gapPos[ i ]
  )

# xb1aDists = xb1aDists[ order( geneOrder ), ]

rowAnnots = data.frame( sapply( genomeNames, function(x)
{
  if ( x %in% speciesClusts[[1]] ) return( "Bx")
  return( "Bo" )
}))
colnames( rowAnnots ) = "species"
colList = list( "species" = c( "Bx" = "red4", "Bo" = "navyblue") )
for ( i in which( numGenes > 1 ) )
{
  outFile = paste0( 
    "results/2020_08_04_consecutive_gene_heat_maps/xb1a_phylo_dists_gene_order_", 
    gapPos[i], '_', gapPos[ i + 1 ], ".png"
    )
  pheatmap::pheatmap(
    xb1aDists[ gapPos[i]:gapPos[ i + 1 ], ],
    cluster_rows = FALSE,
    filename = outFile,
    annotation_col = rowAnnots,
    annotation_colors = colList
    ) 
}

pheatmap::pheatmap(
  sortedMat,
  color = c( bxCols, boCols ),
  cluster_rows = FALSE, 
  gaps_row = gapPos
  )

xb1aDists = t(sapply( 1:length( phyloMats ),
  function(i) phyloMats[[ i ]][ , 16 ]
  ))

pheatmap::pheatmap(
  xb1aDists,
  # color = c( bxCols, boCols ),
  cluster_rows = FALSE,
  annotation
  )

# ------------------------------------------------------------------------------

numBx = length( speciesClusts[[1]] )
numBo = length( speciesClusts[[2]] )
numBxInBxClades = sapply( 1:nrow(sortedMat), 
  function(i) sum( seq( numBx ) %in% sortedMat[ i, 1:numBx ] )
  )

barplot(
  numBxInBxClades[ order(numBxInBxClades) ], 
  col = "red4", 
  border = "red4"
  )
sum( numBxInBxClades == numBx ) / nrow( sortedMat )

# ------------------------------------------------------------------------------

# Make sure I remembered to sort the dist mats 
isSorted = sapply( 2:length(phyloMats), 
  function(i) identical( dimnames( phyloMats[[1]] ), dimnames( phyloMats[[i]]) )
  )
FALSE %in% isSorted

# ---- Calculate the medain distance between alleles from both species ---------

genomeNames = rownames( phyloMats[[1]] )
nGenomes = nrow( phyloMats[[1]] )

CalcMedianDists = function( distMat, isMySpecies )
{
  medianDistVals = sapply( 1:nrow(distMat), function( gIdx )
  {
    if ( isMySpecies[ gIdx ]  ) isMySpecies[ gIdx ] = FALSE
    return( median( distMat[ gIdx, isMySpecies ] ) )
  })
  return( medianDistVals )
}

# Calculate the median distance from each species to every other 
# species
isBx = genomeNames %in% speciesClusts[[1]]
medianBxVals = sapply( 1:length( phyloMats ),
  function(i) CalcMedianDists( phyloMats[[i]], isBx )
  )
medianBoVals = sapply( 1:length( phyloMats ),
  function(i) CalcMedianDists( phyloMats[[i]], !isBx )
  )
bxVals = t( medianBxVals )
boVals = t( medianBoVals )


# ------------------------------------------------------------------------------

bxToBo = sapply( 1:nrow(bxVals), function(i) 
{
  isHgt = sapply( which( isBx ), 
    function(j) bxVals[ i, j ] >= 0.1 && boVals[ i, j ] < 0.1
    )
  return( sum( isHgt ) )
})

boToBx = sapply( 1:nrow(bxVals), function(i) 
{
  isHgt = sapply( which( !isBx ), 
    function(j) bxVals[ i, j ] < 0.1 && boVals[ i, j ] >= 0.1
    )
  return( sum( isHgt ) )
})
sum( bxToBo  )
sum( boToBx )

x = sapply( 1:nrow( bxVals ), 
  function(i) max(  bxVals[ i, ] ) - min( bxVals[ i,  ] )
  )
idx = which.max( x )
idx = which( x == x[ order(x) ][ length(x) / 2 ] )


PlotGeneTree = function( idx )
{
  tree = geneTrees[[ idx ]]
  plot( tree )
  
  speciesCols = sapply( tree$tip.label, function(x) 
  {
    if ( x %in% speciesClusts[[1]] ) return( "red4" )
    return( "navyblue" )
  })
  tiplabels( col = speciesCols, pch = 15 )
}

plotCols = sapply( isBx, function(x) 
{
  if ( x ) return( "red4" )
  return( "navyblue" )
})

# ---- Plot the distance values ------------------------------------------------

plot( 
  medianBxVals,
  medianBoVals,
  col  = plotCols,
  cex  = 1,
  cex.lab = 1.25,
  cex.axis = 1.25,
  pch  = 19,
  xlab = "Median distance to B. xylanisolvens alleles",
  ylab = "Median distance to B. ovatus alleles",
  )
offSetVal  = 0.005
lineWidVal = 7.5
rect(
  xleft   = 0.1 - offSetVal,
  xright  = 0.475 + offSetVal,
  ytop    = 0.1 + offSetVal,
  ybottom = 0.0 - offSetVal,
  lwd     = lineWidVal + 1,
  border  = "gray"
  )
rect(
  xleft   = 0.0 - offSetVal,
  xright  = 0.1 + offSetVal,
  ytop    = 0.475 + offSetVal,
  ybottom = 0.1 - offSetVal,
  lwd     = lineWidVal + 1,
  border  = "gray"
  )
rect(
  xleft   = 0.1 - offSetVal,
  xright  = 0.475 + offSetVal,
  ytop    = 0.1 + offSetVal,
  ybottom = 0.0 - offSetVal,
  lwd     = lineWidVal,
  border  = "red4"
  )
rect(
  xleft   = 0.0 - offSetVal,
  xright  = 0.1 + offSetVal,
  ytop    = 0.475 + offSetVal,
  ybottom = 0.1 - offSetVal,
  lwd     = lineWidVal,
  border  = "darkblue"
  )
legend(
  "topleft",
  legend = c( "B. xylanisolvens allele", "B. ovatus allele" ),
  col    = c( "red4", "navyblue" ),
  pch    = 15
  )

# ------------------------------------------------------------------------------

diffMat = sapply( 1:nrow(medianBxVals), 
  function(i) sapply( 1:ncol( medianBxVals), 
    function(j)  medianBoVals[ i, j ] - medianBxVals[ i, j ]
    )
  )


diffMat[ , !isBx ] = diffMat[ , !isBx ] * -1

pheatmap::pheatmap( diffMat )

PHYLO_DIST_THRESH = 0.1

GetNumSigDiffs = function( inVec )
{
  isSig = abs( inVec ) > PHYLO_DIST_THRESH
  if ( TRUE %in% isSig ) return( sum( inVec[ isSig ] < 0 ) )
  return( 0 )
}

numSigDiffAlleles = sapply( 1:nrow(diffMat),  function(i) 
{
  numSigBx = GetNumSigDiffs( diffMat[ i, isBx ] )
  numSigBo = GetNumSigDiffs( diffMat[ i, !isBx ] )
  return( c( numSigBx, numSigBo ) )
})

isSigDifferent = colSums( numSigDiffAlleles ) != 0

for ( i in  which( isSigDifferent ) ) PlotGeneTree( i )      

pheatmap::pheatmap(
  sortedMat[ isSigDifferent, ],
  color = c( bxCols, boCols )
  )


par( mfrow = c( 2, 1 ), mar = c( 1.1, 4.1, 4.1, 2.1 ) )
barplot( numSigDiffAlleles[1, ], col = "red4", border = "red4" )
barplot( numSigDiffAlleles[2, ], col = "navyblue", border = "navyblue" )

# For bx, the median distance to a bx allele should be smaller than the
# median distance to a bo allele, so find the number of instances where
# the median distance is negative 
numBxMoreLikeBo = sum( diffMat[ , isBx ] < 0 )

# For bo, the opposite
numBoMoreLikeBx = sum( diffMat[ , !isBx ] < 0 )

meanDiffs = rowSums( diffMat ) / nrow( diffMat )
summary( meanDiffs )

numCloserMatch = c( 
  sum( bxToBo > 0 ),
  sum( boToBx > 0)
  )

names( numCloserMatch ) = c(
  "Number of Bx alleles more\nclosely matching Bo alleles",
  "Number of Bo alleles more\nclosely matching Bx alleles"
  )

par( mfrow = c( 1, 1 ), mar = c( 15.1, 4.1, 4.1, 2.1 ) )
barplot(
  numCloserMatch,
  col = c("red4", "navyblue"),
  border = c("red4", "navyblue"),
  pch = 19,
  # xlab = "Median distance to Bx alleles",
  ylab = "Number of alleles with closer match to the other species",
  las = 3,
  ylim = c( 0, 50)
  )
legend(
  "topright",
  legend = c("Bx", "Bo"),
  col    = c("red4", "navyblue"),
  pch    = 15
  )

# ------------------------------------------------------------------------------

save( file = rData, list = ls() )

# ------------------------------------------------------------------------------