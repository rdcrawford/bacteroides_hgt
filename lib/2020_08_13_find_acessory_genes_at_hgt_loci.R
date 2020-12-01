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

load("data/2020_08_13_find_hgt_accessory_genes.Rdata")

# Make a list with the cluster indexes of each gene  
clIdxs = sapply( 1:length(genePtr$clustList), 
  function(i) list( rep( i , length( genePtr$clustList[[i]] ) ) )
  )

# Create a data-frame with the genenames, genome ids, and cluster indexes for
# easy look up
coreGeneData = cbind.data.frame(
  unlist( genePtr$clustList ),
  sapply(unlist( genePtr$clustList ), GetGenomeId ),
  unlist( clIdxs ),
  stringsAsFactors = FALSE
  )

# Create a data-frame with the meta data on the genes 
geneData = GetGeneData()

# For each core gene, look up the name of the cluster in the entire cluster
# list
coreGeneClIdxs = GetWholeGenomeClust()

# Get the number of genomes
nGenomes = length( genomeNames )

# Get the genome names of each gff file
gfNames = sapply( 1:length( genePtr$gfList ),
  function(i)  GetGenomeId( genePtr$gfList[[i]][ 1, 1 ] )
  )

# Update the gff files to include the cluster name for each gene
gfList = sapply( 1:length( genePtr$gfList ),
  function(i) list( AddClIdxsToGff( genePtr$gfList[[i]], gfNames[i] ) )
  )

# Sort the gff files to match the genome names
genomeOrder = sapply( genomeNames, function(x) which( gfNames == x ) )
gfList = gfList[ order( genomeOrder ) ]

# ------------------------------------------------------------------------------

# Remove redundant genomes from the analysis
dupGenomes = c(
  "Bacteroides_xylanisolvens_CL03T12C04",
  "Bovatus_CL03T12C18",
  "Bxylanisolvens_CL03T12C04"
  )
isKeeper = !genomeNames %in% dupGenomes
bxVals   = bxVals[ , isKeeper ] 
boVals   = boVals[ , isKeeper ]
isBx     = isBx[ isKeeper ]
gfList   = gfList[ isKeeper ]
nGenomes = sum( isKeeper )

# ------------------------------------------------------------------------------

# Find genes share at hgt loci in bx
bxHgtGenes = FindAccGenesAtHgtLoci( bxVals, boVals, isBx )

# Find genes at hgt loci in bo
boHgtGenes = FindAccGenesAtHgtLoci( boVals, bxVals, !isBx )

# ------------------------------------------------------------------------------

for ( i in 1:length(bxHgtGenes) )
{
  rGenome = GetRefGenomeIdx( names( bxHgtGenes )[i] )
  if ( is.matrix( bxHgtGenes[[i]] ) ) 
    MakeLocusGeneHeatmap( bxHgtGenes[[i]], rGenome )
}
i = 4
MakeLocusGeneHeatmap( bxHgtGenes[[i]], rGenome )
# ------------------------------------------------------------------------------

# Get the co-phylogenetic distances for the core gene tree
coPhyloDists = GetCoreGenePyloMat()

#


# ------------------------------------------------------------------------------

save( file = "data/2020_08_13_find_hgt_accessory_genes.Rdata", list = ls())
save( 
  file = "~/Desktop/2020_08_13_find_hgt_accessory_genes.Rdata", 
  list = ls()
  )
# ------------------------------------------------------------------------------