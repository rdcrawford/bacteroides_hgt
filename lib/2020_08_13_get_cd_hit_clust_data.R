# ------------------------------------------------------------------------------
# 
# 2020/08/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( Rcpp )
library( cognac )
library( ape )
library( maps )
library( phytools )
library( RColorBrewer )

# ------------------------------------------------------------------------------

HGT_THRESH    = 0.1
G_START       = 4
G_END         = 5
G_CLUST       = 7
MAX_GENE_DIST = 5000

# ------------------------------------------------------------------------------

load( "2020_07_30_analyze_tree_order.Rdata" )

# ------------------------------------------------------------------------------

AddClustIdxToGff = function( geneEnv )
{
  for ( j in 1:length(geneEnv$gfList ) )
  {
    clIdxs = vector( "integer", nrow(geneEnv$gfList[[j]]) )

    for ( i in 1:length(geneEnv$clustList) )
    {
      isClustGene =
        geneEnv$gfList[[j]][ , FEATURE_ID ] %in% geneEnv$clustList[[i]]

      if ( TRUE %in% isClustGene )
      {
        for ( k in which(isClustGene) ) clIdxs[ k ] = i
      }
    }

    geneEnv$gfList[[j]] = cbind.data.frame(
      geneEnv$gfList[[j]],
      clIdxs,
      stringsAsFactors = FALSE
      )
  }
}

GetCoreGeneClustIdxs = function( geneEnv )
{
  
  cgClustList = strsplit( algnData$geneData$clGeneIds, ',' )
  sapply( cgClustList, function(x) which( geneEnv$clustList == x ) )
}


GetWholeGenomeClust = function()
{
  cgClList = strsplit( algnData$geneData$clGeneIds, ',' )
  
  nGenes = length( geneEnv$clustList )
  cgClIdxs = sapply( 1:length( cgClList ), function(i)
  {
    j = 1
    while( 
      FALSE %in% ( cgClList[[i]] %in% geneEnv$clustList[[j]] ) && 
      j <= nGenes 
      )
    {
      j = j + 1
    }
    if ( j < nGenes ) return( j )
    return( NA )
  })
  
  return( cgClIdxs )
}

# ------------------------------------------------------------------------------

# Parse the cd-hit data
cdHitClstFile = paste0(
  "analysis/2019_12_06_core_gene_tree/temp_congnac_files/",
  "cdHitClusters.faa.clstr"
  )
cognac::ParseCdHit(
  cdHitClstFile = cdHitClstFile,
  isBinary      = TRUE,
  clSizeThesh   = 1,
  geneEnv       = geneEnv
  )

# Make a list with the cluster indexes of each gene  
clIdxs = sapply( 1:length(geneEnv$clustList), 
  function(i) list( rep( i , length( geneEnv$clustList[[i]] ) ) )
  )

# For each core gene, look up the name of the cluster in the entire cluster
# list
coreGeneClIdxs = GetWholeGenomeClust()

# Get the number of genomes
nGenomes = length( genomeNames )

# Get the genome names of each gff file
gfNames = sapply( 1:length( geneEnv$gfList ),
  function(i)  GetGenomeId( geneEnv$gfList[[i]][ 1, 1 ] )
  )

# Update the gff files to include the cluster name for each gene
gfList = sapply( 1:length( geneEnv$gfList ),
  function(i) list( AddClIdxsToGff( geneEnv$gfList[[i]], gfNames[i] ) )
  )

# Sort the gff files to match the genome names
genomeOrder = sapply( genomeNames, function(x) which( gfNames == x ) )
gfList = gfList[ order( genomeOrder ) ]

geneMat = t( geneEnv$geneMat )

# ------------------------------------------------------------------------------

save( file = "data/2020_08_13_get_cd_hit_clust_data.Rdata", list = ls())

# ------------------------------------------------------------------------------