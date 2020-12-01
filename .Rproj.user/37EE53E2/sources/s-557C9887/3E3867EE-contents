# ------------------------------------------------------------------------------
# Run cognac
# 2020/07/27
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This script creates a concatenated gene alignment with cognac. The alignment
# is then split into the component genes and 
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( cognac )
library( ape )
library( phytools )

# ------------------------------------------------------------------------------

# Split an alignment into paritions
Rcpp::sourceCpp("lib/R/GetSeqPatitions.cpp")

# Divide the tree into n clusters based off of distance
source( "lib/R/GetClustList.R" ) 

# ------------------------------------------------------------------------------

fastaFiles = system( "ls data/fasta/*.fasta", intern = TRUE )
gffFiles   = system( "ls data/gene_annotations/*.gff", intern = TRUE )
genomeIds  = sapply( fastaFiles, ExtractGenomeNameFromPath )

# Remove redundant genomes from the analysis
dupGenomes = c(
  "Bacteroides_xylanisolvens_CL03T12C04",
  "Bovatus_CL03T12C18",
  "Bacteroides_ovatus_CL02T12C04"
  )
isKeeper   = !genomeIds %in% dupGenomes
fastaFiles = fastaFiles[ isKeeper ]
gffFiles   = gffFiles[ isKeeper ]
genomeIds  = genomeIds[ isKeeper ]

outDir = "analysis/2019_12_06_core_gene_tree/"

rData = "data/2020_09_23_cognac_run_data.Rdata"

# ---- Create the core gene tree -----------------------------------------------

geneEnv = CreateGeneDataEnv( gffFiles, fastaFiles, genomeIds, outDir )

algnData = cognac(
  fastaFiles    = fastaFiles,
  featureFiles  = gffFiles,
  revTranslate  = TRUE,
  threadVal     = 1,
  outDir        = outDir,
  keepTempFiles = TRUE
  )
save( file = rData, list = ls() )

# ---- Get the individual gene alinments ---------------------------------------

# Read in the fasta file with the amino acid alignment
algn = ParseFasta( algnData$aaAlgnPath )


# Make a numeric vector with the end of the genes in the alignment
aaPartitions = sapply( algnData$geneData$aaGenePartitions, 
  function(x) as.numeric( strsplit( x, '-' )[[ 1 ]][ 2 ] )
  )

# Split the alignment into the core genes. Store the character
# vector containing each gene alignment as an element in a list
algnList = GetSeqPatitions( algn, aaPartitions )

# ---- Run fastTree -----------------------------------------------------------

# Create vectors with the paths to the alignments and trees 
algnPaths = sapply( seq( length( algnList ) ),
  function(i) paste0(outDir, "aa_gene_algns/", "gene_", i, "_algn.fasta")
  )
treePaths = sapply( algnPaths, function(x) gsub( "fasta", "tre", x ) )

# For each alignment...
for ( i in seq( length(algnList ) ) )
{
  # Write the alignment to a fasta file
  sink( algnPaths[i] )
  for ( j in 1:length( algnList[[i]] ) )
    cat('>', genomeIds[j], '\n', algnList[[i]][j],'\n', sep = '')
  sink()

  # Run fastTree
  system( paste( "fasttree < ", algnPaths[i], '>', treePaths[i] ) )
}

# Read in the gene trees
geneTrees = sapply( treePaths, 
  function(x) list( midpoint.root( read.tree( x ) ) )
  )

# Get the phylogenetic distances
phyloMats = sapply( 1:length(geneTrees),
  function(i) list( cophenetic.phylo( geneTrees[[i]] ) )
  )
save( file = rData, list = ls() )

# ------------------------------------------------------------------------------

# Run fasttree to get an ~ML tree for both the amino acid alignment and 
# the reverse translated nuclueotide alignment 
cognacTree = paste0( outDir, "cognac.tre" )
system( paste0( "fasttree -nt < ", algnData$ntAlgnPath, " > ", cognacTree ) )
cognacAaTree = paste0( outDir, "cognac_aa.tre" )
system( paste0( "fasttree < ", algnData$aaAlgnPath, " > ", cognacAaTree ) )

# Plot the trees 
cognacTree = read.tree( cognacTree )
cognacTree = midpoint.root( cognacTree )
plot( cognacTree )

cognacAaTree = read.tree( cognacAaTree )
cognacAaTree = midpoint.root( cognacAaTree )
plot( cognacAaTree )

# ------------------------------------------------------------------------------

# Find the species of each genome
speciesClusts = GetClustList( cognacAaTree, 2 ) 

GetSpeciesLabels = function( genomeIds, speciesClusts )
{
  species = sapply( genomeIds, function(x)
  {
    if ( x %in%  speciesClusts[[1]] ) return( "B. xylanisolvens" )
    return( "B. ovatus" )
  })
  return( species )
}


# ---- Save the data -----------------------------------------------------------

save( file = rData, list = ls() )

# ------------------------------------------------------------------------------