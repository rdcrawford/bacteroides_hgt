# ------------------------------------------------------------------------------
# Find Accessory Genes At HGT Loci Functions
# 2020/08/19
# Ryan D. Crawford
# ------------------------------------------------------------------------------

# Find genes exhibiting patterns of HGT in the data-set
FindHgtGenes = function( refVals, qryVals, isMySpecies )
{
  numGenomes  = sum( isMySpecies )
  numGenes    = nrow( refVals )
  speciesIdxs = which( isMySpecies )
  isHgt       = matrix( FALSE, nrow = numGenes, ncol = numGenomes )

  for ( i in 1:nrow( refVals ) )
  {
    for ( j in seq( numGenomes ) )
    {
      cIdx = speciesIdxs[j]
      if ( refVals[ i, cIdx ] > 0.1 && qryVals[ i, cIdx ] <= 0.1 )
        isHgt[ i, j ] = TRUE
    }
  }
    
  return( isHgt )
}

# Find the core gene names of genes with patterns of HGT
GetHgtCgClstIdxs = function( isHgt )
{
  isHgtGene  = sapply( 1:nrow( isHgt ), function(i) TRUE %in% isHgt[ i, ] )
  return( which( isHgtGene ) )
}

# Find the genomes that have HGT genes
GetHgtGenomes = function( isMySpecies, isHgt, cgHgtClusts )
{
  genomeIdxs = which( isMySpecies )
  return( sapply( cgHgtClusts, function(i) list( genomeIdxs[ isHgt[i, ] ] ) ) )
}

# increment the index of the current gene in the gff file
GetNextGene = function( inVal )
{
  return( inVal + 1 )
}

# decrement the index of the current gene in the gff file
GetPrevGene = function( inVal )
{
  return( inVal - 1 )
}

FindStrandDir = function( r, q, rIdx, qIdx )
{
  # Find the number of genes in the same direction as the reference
  # and the number of genes in the opposite direction of the reference in 
  # the query
  stepVal = 0
  while ( 
    gfList[[ rIdx ]][ r - stepVal, G_CLUST ] == 
    gfList[[ rIdx ]][ r + stepVal, G_CLUST ]
    )
  {
    stepVal = stepVal + 1
  }
  
  # 
  lhsRefGene = gfList[[ rIdx ]][ r - stepVal, G_CLUST ]
  rhsRefGene = gfList[[ rIdx ]][ r + stepVal, G_CLUST ]
  lhsQryGene = gfList[[ qIdx ]][ q - stepVal, G_CLUST ]
  rhsQryGene = gfList[[ qIdx ]][ q + stepVal, G_CLUST ]

  # Count the number of genes in the same direction and the opposite direction
  isSameDir = ( lhsRefGene == lhsQryGene ) || ( rhsRefGene == rhsQryGene )
  isOppDir  = ( lhsRefGene == rhsQryGene ) || ( rhsRefGene == lhsQryGene )
  
  # If there are no matches on either side, return null indicating the 
  # strand direction cannot be assigened
  if ( !isSameDir && !isOppDir ) return( NULL )
  
  # if thereare more genes in the opposite direction return false, true if
  # they are in the same direction
  if ( isSameDir ) return( TRUE )
  return( FALSE )
}

# If the genes are in the same direction, make the function to
# # get the next gene the same
GetSharedContigGenes = function( rIdx, r, qIdx, q, GetRefGeneFn, GetQryGeneFn )
{
  while (
    gfList[[ rIdx ]][ GetRefGeneFn( r ), G_CLUST ] == 
    gfList[[ qIdx ]][ GetQryGeneFn( q ), G_CLUST ]
    )
  {
    r = GetRefGeneFn( r )
    q = GetQryGeneFn( q )
  }

  return( r )
}

# For each index given find the indicies in the regference genome (as 
# indicated by gIdx). Return a matrix with the indicies in the gff file
# indicating the left and right most positions shared between each  
# pair of genomes.
FindSyntGenesAtLocus = function( gIdx, toTest, clustIdxs )
{
  # Get the row index of the gene in the reference genome
  r = clustIdxs[ gIdx ]
  
  # Iterate over the members of the other species and ding the maximium
  # number of genes shared between the genome of interest and the members
  # of the opposite species
  syntGeneIdxs = sapply( toTest, function(i)
  {
    # Get the row index in the gff file of the gene in the reference genome 
    q = clustIdxs[ i ]

    # Find if the genes are in the direction as the reference genome
    isSameDir = FindStrandDir( r, q, gIdx, i )
    
    # If there are no shared genes on either side, return 0
    if ( is.null( isSameDir ) )
    {
      return( c( clustIdxs[ gIdx ], clustIdxs[ gIdx ] ) )
    }
    else if ( isSameDir )
    {
      leftIdx  = GetSharedContigGenes( gIdx, r, i, q, GetPrevGene, GetPrevGene )
      rightIdx = GetSharedContigGenes( gIdx, r, i, q, GetNextGene, GetNextGene )
    }
    else 
    {
      leftIdx  = GetSharedContigGenes( gIdx, r, i, q, GetPrevGene, GetNextGene )
      rightIdx = GetSharedContigGenes( gIdx, r, i, q, GetNextGene, GetPrevGene )
    }
    return(  c( leftIdx, rightIdx ) )
  })

  return( syntGeneIdxs )
}

# Find the maximium number of genese shared between the 
# referece genome (indicated by gene index) and the genomes of the 
# opposite species. Once the maximium numer of shared genes is identified, 
# return a list with: the genome, the gene used as the seed, and the 
# names of the genes which were transferred
FindSyntenicHgtLoci = function( clustIdxs, gIdx, isMySpecies )
{
  # Get the indexes of the opposite species
  toTest = which( !isMySpecies )

  # Find the farmost right and left gene where there is still a match
  # in the opposite species
  syntGeneIdxs = FindSyntGenesAtLocus( gIdx, toTest, clustIdxs )
  
  # Calculate the number of genes identified in each HGT event
  nGenes = sapply( 1:ncol( syntGeneIdxs ), 
    function(j) syntGeneIdxs[ 2, j ] - syntGeneIdxs[ 1, j ] + 1
    )
  
  # Get the maximium number of genes shared between this genome and
  # a genome of the opposite species
  maxGenes = which.max( nGenes )
  startPos = syntGeneIdxs[ 1, maxGenes ]
  endPos   = syntGeneIdxs[ 2, maxGenes ]
  hgtGenes = gfList[[  gIdx  ]][ startPos:endPos, G_CLUST ]
  
  # Look up the gene name of the gene used as the seed
  seedGene = gfList[[ gIdx ]][ clustIdxs[ gIdx ], G_CLUST ]

  # Create a list to output where the first evlement is the genome 
  # name and the second element is the syntenic co-linear genes
  # shared between this genome and the opposite species
  return( list( gIdx, seedGene, hgtGenes ) )
}

# Iterate over the hgt events and identify any instances where 
# the exact set of genes was identified in two genomes. If any cases
# of this were identified, merge the genome IDs for the event. The 
# de-duplicated set of HGT events is returned.
RemoveDuplicateHgtEvents = function( hgtEventData )
{
  # Look up all of the gene at the same locus as the hgt gene
  locGenes = sapply( 1:length( hgtEventData ), function(j)
  {
    genes = unique( as.vector( hgtEventData[[ j ]][[ HGT_GENES ]] ) )
    return( list( genes[ order( genes ) ] ) )
  })
  
  # Idenfify any duplicated gene sets 
  isDuplicate = duplicated( locGenes )

  # If there are no duplicates, there's nothing to do return the 
  # list unchanged
  if ( !TRUE %in% isDuplicate ) return( hgtEventData )
  
  # Create the names for the output list
  for ( i in which(!isDuplicate) )
  {
    # Find the identical gene sets to the index
    isMatch = sapply( 1:length( locGenes ), function(j) 
      identical( locGenes[[j]], locGenes[[i]] )
      )
    
    # Find the genome names and the seed gene names and 
    # collapse them for the 
    hgtGenomeIds = unique( sapply( which( isMatch ), 
      function(i) hgtEventData[[i]][[ GENOME_NAME ]]
      ))
    hgtSeedGenes = unique( sapply( which( isMatch ), 
      function(i) hgtEventData[[i]][[ CLUST_IDX ]]
      ))
    hgtEventData[[i]][[ CLUST_IDX ]]   = hgtSeedGenes
    hgtEventData[[i]][[ GENOME_NAME ]] = hgtGenomeIds
  }
  return( hgtEventData[ !isDuplicate ] )
}

# This functon iterates over the hgt event data and finds any instances
# where the same genes were identified within another set of genes, allowing
# for one mismatch. Nested HGT events are merged. The hgt events list, 
# sorted from largest to smallest is returned
FindNestedHgtEvents = function( hgtEventData )
{
  # Find the number of hgt genes in sample
  nEvents    = length( hgtEventData )
  eventRange = seq( nEvents )
  numGenes   = sapply( eventRange, 
    function(i) length( hgtEventData[[i]][[ HGT_GENES ]] )
    )
  hgtEventData = hgtEventData[ order( numGenes, decreasing = TRUE ) ]
  
  # Create a logical vector to indicate whether an hgt event is still
  # in the analysis
  isClassified = vector( "logical", nEvents )
  
  # Initialize the reference and querry 
  rIdx = 1
  qIdx = 2
  toKeep = rIdx
  isClassified[ rIdx ] = TRUE
  
  # While there are events that remain unclassified...
  while ( FALSE %in% isClassified )
  {
    # Get the indexes of the unclassified evnets
    isUnclassified = which( !isClassified )
    
    # Find if all of the genes in the querry are present in the reference
    isInRef = hgtEventData[[ qIdx ]][[ HGT_GENES ]] %in% 
      hgtEventData[[ rIdx ]][[ HGT_GENES ]]
    
    # Find if the HGT events share a marker gene
    hasSharedMarker = hgtEventData[[ qIdx ]][[ CLUST_IDX ]] %in% 
      hgtEventData[[ rIdx ]][[ CLUST_IDX ]]
    
    # Get a count of the number of missing genes 
    numMissing = sum( !isInRef )
    
    # If all except one of these genes are present in the reference, 
    # merge the events and they have a shared phylogenetic marker with
    # evidence of HGT
    if ( numMissing <= 1 && ( TRUE %in% hasSharedMarker ) )
    {
      # Merge the genome ids
      hgtEventData[[ rIdx ]][[ GENOME_NAME ]] = unique(c(
        hgtEventData[[ rIdx ]][[ GENOME_NAME ]], 
        hgtEventData[[ qIdx ]][[ GENOME_NAME ]]
        ))
  
      # Add the missing gene
      if ( numMissing == 1 )
      {
        # Get the index of where the gene will be added
        missingIdx = which( !isInRef ) - 1
        
        # Find the gene to add
        toAdd = hgtEventData[[ qIdx ]][[ HGT_GENES ]][ !isInRef ]
        
        # Append the new gene to the vector
        hgtEventData[[ rIdx ]][[ HGT_GENES ]] = append(
          hgtEventData[[ rIdx ]][[ HGT_GENES ]], 
          toAdd,
          missingIdx
          )
      }
      
      # Mark this query as classified
      isClassified[ qIdx ] = TRUE
    }
    
    # If this is the last query set a new references
    if ( qIdx >= isUnclassified[ length( isUnclassified ) ] )
    {
      # Mark this event as classified
      rIdx = isUnclassified[ 1 ]
      isClassified[ rIdx ] = TRUE
      toKeep[ length( toKeep ) + 1 ] = rIdx
      
      # Reset the query to the first unclassified sequence after the reference
      if ( length( isUnclassified ) > 2 ) qIdx = isUnclassified[ 2 ]
      
    # Set a new query
    } else {
      
      qIdx = isUnclassified[ isUnclassified > qIdx ][ 1 ]
    }
  }

  return( hgtEventData[ toKeep ] )
}

# This function initializes the list of hgt event data. For each hgt gene,
# and for each HGT genome with the putative HGT gene, a list created by 
# the "FindSyntenicHgtLoci" function is created. Each list element is composed
# of a three part list contining: the genome id, the name of the 
# gene used as the seed, and the names of the genes that are perfectly
# syntenic and colinear between this genome and at least one other genome
# of the opposite species. 
FindGenesAtHgtLoci = function( hgtClusts, hgtGenomes, isMySpecies )
{
  # Create an empty list to store the
  hgtEventData = vector( "list", 0 )

  # For each of the HGT genes...
  for ( clustIdx in 1:length( hgtClusts ) )
  {
    # Find the position in the gff file corresponding to the hgt
    # cd-hit cluster name in each genome
    clustIdxs = sapply( 1:length( gfList ),
      function(i) which( gfList[[i]][ , G_CLUST ] == hgtClusts[ clustIdx ] )
      )

    # A list of the syntenic and colinear genes shared between each 
    # genome with a putative HGT event
    locusHgtList = sapply( hgtGenomes[[ clustIdx ]],
      function(j) list( FindSyntenicHgtLoci( clustIdxs, j, isMySpecies ) )
      )

    # Find the number of genes with signitures of hgt
    numGenes = length( locusHgtList )

    # Append the results to the list of matricies to output
    if ( numGenes > 1 )
    {
      listRange =
        ( length( hgtEventData ) + 1 ):( length( hgtEventData ) + numGenes )
    } else {
      listRange = length( hgtEventData ) + 1
    }
    hgtEventData[ listRange ] = locusHgtList
    
  }
  return(  hgtEventData )
}

# Calcuate the Pr HGT gene | no HGT evidence
GetHgtGeneProbs = function( hgtEventData, isMySpecies, specNames )
{
  # Find the species that are not
  spIdxs  = which( isMySpecies )
  otherSp = which( !isMySpecies )
  
  for ( i in 1:length( hgtEventData ) )
  {
    if ( length( hgtEventData[[i]][[ ACC_GENES ]] ) )
    {
      genes      = hgtEventData[[i]][[ ACC_GENES ]]
      nonHgtIdxs = spIdxs[ !spIdxs %in% hgtEventData[[i]][[ GENOME_NAME ]] ]
         
      prMat = rbind( 
        CalcGenePrs( genes, nonHgtIdxs ), 
        CalcGenePrs( genes, otherSp )
        )
      row.names( prMat ) = specNames
      hgtEventData[[i]][[ ACC_GENE_PRS ]] = prMat
    }
    else
    {
      hgtEventData[[i]][[ ACC_GENE_PRS ]] = numeric()
    }
  }
  
  return( hgtEventData )
}

# For each of the genes calcuate the probability of each gene for the input 
# genomes
CalcGenePrs = function( genes, genomes )
{
  myGenes = geneMat[ genes, genomes ]
  if ( is.null( dim( myGenes ) ) ) return( sum( myGenes ) / length( genomes ) )
  return( rowSums( myGenes ) / length( genomes ) )
}

# For a genome, look up the annotations for a gene from the gff 
# file
GetGeneAnnotations = function( geneNames, genomeIdx )
{
  rowIdxs = sapply( geneNames, function(x) 
  {
    isMyGene = gfList[[ genomeIdx ]][ , G_CLUST ] == x
    return( which( isMyGene )[1] )
  })
  return( gfList[[ genomeIdx ]][ rowIdxs, G_ANNOT ] )
}

# Append the gene annotations to each list describing the 
# hgt event
LookUpHgtGeneAnnots = function( hgtEventData )
{
  for ( i in 1:length( hgtEventData ) )
  {
    hgtEventData[[i]][[ GENE_ANNOTS ]] = GetGeneAnnotations(
      hgtEventData[[i]][[ HGT_GENES ]],
      hgtEventData[[i]][[ GENOME_NAME ]][1]
      )
  }
  return( hgtEventData )
}

# Print the information on a single HGT event
WriteHgtData = function( hgtEvent )
{
  # Find if any of the genes contribute to pan genome expansion
  hgtGenes   = hgtEvent[[ HGT_ACC_GENES ]]
  nAccGenes  = length( hgtGenes )
  nGenomes   = length( hgtEvent[[ GENOME_NAME ]] )
  nGenes     = length( hgtEvent[[ HGT_GENES ]] )
  nCoreGenes = length( hgtEvent[[ CLUST_IDX ]])
  
  # Write the names of the genomes where the event occured
  cat("  -- Present in ", nGenomes, " Genomes:\n", sep = '')
  for ( i in hgtEvent[[ GENOME_NAME ]] )
  {
    cat("     -- ", genomeNames[i], '\n', sep = '')
  }
  
  # Get the 
  cat("  -- ", nGenes, " syntenic, colinear genes with ", nCoreGenes,
    " high probability HGT core genes and ", 
    nAccGenes, " high probability HGT accessory genes","\n", sep = ''
    )
  for ( i in 1:length( hgtEvent[[ HGT_GENES ]] ) )
  {
    gene = hgtEvent[[ HGT_GENES ]][i]
    isCg = colnames( hgtEvent[[ CG_DISTS ]] ) == gene
    cat( "     -- ", gene, ": ", hgtEvent[[ GENE_ANNOTS ]][i], sep = '' )
    if ( TRUE %in% isCg )
    {
      if ( gene %in% hgtEvent[[ CLUST_IDX ]] )
      {
        cat( "; HGT core gene" )
      } 
      else 
      {
        cat( "; Core gene" )
      }
      
      cat( "; " )
      for ( i in 1:nrow( hgtEvent[[ CG_DISTS ]] ) )
      {
        cat( 
          row.names( hgtEvent[[ CG_DISTS ]] )[i], ": ", 
          hgtEvent[[ CG_DISTS ]][ i, isCg ],
          sep = ''
          )
        if ( i != nrow( hgtEvent[[ CG_DISTS ]] ) ) cat(", ")
      }
    }
    else if ( gene %in% hgtGenes )
    {
      cat( "; High confidence HGT accessoy gene")
    } 
    
    cat( "\n" )
  }
}

# For each HGT event in the list, write the data to the input file path
WriteHgtEvents = function( hgtEventData, outFile )
{
  sink( outFile )
  for ( i in 1:length( hgtEventData ) )
  {
    cat("HGT Event:", i, '\n')
    WriteHgtData( hgtEventData[[i]] )
  }
  sink()
}

# The input matricies are analyzed to identify core genes with 
# evidence of HGT (from query to reference). Reference is the recipient 
# of the HGT sequences and Query is the potential donor of the HGT sequences.
# Once HGT core genes are identified, accessory genes at HGT loci are 
# identified. A list of hgt is output with nested and duplicated events 
# collapsed. 
FindAccGenesAtHgtLoci = function( 
  refVals,     # Geneome x gene matrix with the distances to the ref allels
  qryVals,     # Geneome x gene matrix with the distances to the query allels
  isMySpecies, # Logical vector indicating whether a genome is the ref species
  specNames    # Char vector length 2 with the names of the species (ref first)
  )
{
  # Find the genes and alleles for which there is evidence of HGT
  isHgt = FindHgtGenes( refVals, qryVals, isMySpecies )

  # Get the cluster names of each gen with HGT
  cgHgtClusts = GetHgtCgClstIdxs( isHgt )

  # Convert the name of the core genome clusuters to the whole genome clusers
  hgtClusts = coreGeneClIdxs[ cgHgtClusts ]

  # Find the genomes for which the hgt evenst occurs
  hgtGenomes = GetHgtGenomes( isMySpecies, isHgt, cgHgtClusts )
  
  # Create a list of the syntenic genes present in both species at each 
  # hgt locus
  hgtEventData = FindGenesAtHgtLoci( hgtClusts, hgtGenomes, isMySpecies )
  
  # Remove any duplicates where multiple genes or multiple genomes have the
  # synteny and colineary at each hgt locus
  hgtEventData = RemoveDuplicateHgtEvents( hgtEventData )

  # Find any instances where gene content at the locus differers but 
  hgtEventData = FindNestedHgtEvents( hgtEventData )
  
  # Add the distances to the core genes in this HGT event
  hgtEventData = GetCoreGeneDists( hgtEventData, refVals, qryVals, specNames )
  
  # Find the hgt genes at each locus
  hgtEventData = FindAccHgtGenes( hgtEventData )
  
  # Calculate the probability of accessory genes at HGT loci
  hgtEventData = GetHgtGeneProbs( hgtEventData, isMySpecies, specNames )
  
  # Look up new gene flow into the species
  hgtEventData = FindNewHgtGenes( hgtEventData )
  
  # Look up the gene annotations
  hgtEventData = LookUpHgtGeneAnnots( hgtEventData )
  
  # Find the nearest allele matching the hgt alleles
  return( hgtEventData )
}

# Identify the accesspry genes that are only present in the presence of
# core genome variation characteristic of HGT
FindNewHgtGenes = function( hgtEventData )
{
  for ( i in 1:length( hgtEventData ) )
  {
    if ( length( hgtEventData[[i]][[ ACC_GENE_PRS ]] ) )
    {
      # Find the genes with a zero probability of being present in the
      # same species in the absence of an hgt event
      isNewAccGene = hgtEventData[[i]][[ ACC_GENE_PRS ]][ 1, ] == 0
      hgtAccGenes  = hgtEventData[[i]][[ ACC_GENES ]][ isNewAccGene ]
    }
    else
    {
      hgtAccGenes = numeric()
    }
    
    # Add the genes to the hgt event 
    hgtEventData[[i]][[ HGT_ACC_GENES ]] = hgtAccGenes
  }
  return( hgtEventData )
}

# Get the distances to the core genes involved in each HGT event
GetCoreGeneDists = function( hgtEventData, refVals, qryVals, specNames )
{
  for ( i in 1:length( hgtEventData ) )
  {
    # Find the core genes present in this hgt event using the 
    # first genome as the representitive
    isEventCg = hgtEventData[[i]][[ HGT_GENES ]] %in% coreGeneClIdxs
    eventCgs  = hgtEventData[[i]][[ HGT_GENES ]][ isEventCg ]
    gIdx      = hgtEventData[[i]][[ GENOME_NAME ]][1]
    
    # Create a matrix where the first 
    distMat = sapply( eventCgs, function(x) 
    {
      rIdx = which( coreGeneClIdxs == x )
      return( c( refVals[ rIdx, gIdx ], qryVals[ rIdx, gIdx ] ) )
    })
    colnames( distMat )  = eventCgs
    row.names( distMat ) = specNames
   
    # Append the matrix to the 
    hgtEventData[[i]][[ CG_DISTS ]] = distMat
  }
  return( hgtEventData )
}

# Look up the genes corresponding to accessory gene in the HGT data
FindAccHgtGenes = function( hgtEventData )
{
  
  for( i in 1:length( hgtEventData ) )
  {
    # Find the core genes amound the putative 
    hgtGenes   = hgtEventData[[i]][[ HGT_GENES ]]
    isCoreGene = hgtGenes %in% colnames( hgtEventData[[i]][[ CG_DISTS ]] )
    
    # Add the genes to the list
    hgtEventData[[i]][[ ACC_GENES ]] = hgtGenes[ !isCoreGene ]
  }
  
  return( hgtEventData )
}

# Count the number of HGT accessory genes at each locus
GetNumHgtGenes = function( hgtEventData )
{
  numGenes = sapply( 1:length( hgtEventData ),
    function(i) length( hgtEventData[[i]][[HGT_ACC_GENES]] )
    )
  return( numGenes )
}

# Create a character vector of HGT accessory genes at each locus
GetHgtGenes = function( hgtEventData )
{
  genes = unlist(sapply( 1:length( hgtEventData ),
    function(i) hgtEventData[[i]][[ HGT_GENES ]]
    ))
  return( genes )
}

# Plot a gene by genome heatmap for a given hgt locus
MakeGeneHeatmap = function( hgtEvent )
{
  
  hgtMat              = geneMat[ hgtEvent[[ HGT_GENES ]],  ]
  colnames( hgtMat )  = genomeNames
  row.names( hgtMat ) = hgtEvent[[ GENE_ANNOTS ]]
  
  pheatmap::pheatmap(
    hgtMat[ , c( which(isBx), which(!isBx) ) ],
    color             = c( "grey", "red4"),
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    # gaps_col          = gaps,
    # main              = hmTitle,
    annotation_col    = rowAnnots,
    annotation_colors = colList
    )
}

# Filter an hgt event to the first hgt accessory gene or first hgt core gene
FilterHgtEventData = function( hgtEvent )
{
  # Find the instances where there may be HGT
  isHgt = which(sapply( 1:ncol( hgtEvent[[ CG_DISTS ]] ),
    function(j) 
      hgtEvent[[ CG_DISTS ]][ 1, j ] > hgtEvent[[ CG_DISTS ]][ 2, j ]
    ))
  
  # Find the positions of the first and last gene that has evidence of 
  # HGT
  minHgt    = min( isHgt )
  maxHgt    = max( isHgt )
  startGene = colnames( hgtEvent[[ CG_DISTS ]] )[ minHgt ]
  endGene   = colnames( hgtEvent[[ CG_DISTS ]] )[ maxHgt ]
  
  if ( minHgt == maxHgt )
  {
    hgtEvent[[ CG_DISTS ]] = matrix( hgtEvent[[ CG_DISTS ]][ , maxHgt ] )
    colnames( hgtEvent[[ CG_DISTS ]] ) = startGene
  } else {
    hgtEvent[[ CG_DISTS ]] = hgtEvent[[ CG_DISTS ]][ , minHgt:maxHgt ]
  }
  
  
  # Find the first and last hgt gene at this locus
  startPos = which( hgtEvent[[ HGT_GENES ]] == startGene )
  endPos   = which( hgtEvent[[ HGT_GENES ]] == endGene )
  
  # If there are any hgt accessory genes, find if there are 
  # andy accessory genes that are 
  if ( length( hgtEvent[[ HGT_ACC_GENES ]] ) )
  {
    isHgtAccGene = hgtEvent[[ HGT_GENES ]] %in% hgtEvent[[ HGT_ACC_GENES ]]
    firstAccGene = min( which( isHgtAccGene ) )
    lastAccGene  = max( which( isHgtAccGene ) )
    if ( firstAccGene < startPos ) startPos = firstAccGene
    if ( lastAccGene > endPos ) endPos = lastAccGene
  }
  hgtEvent[[ HGT_GENES ]]   = hgtEvent[[ HGT_GENES ]][ startPos:endPos ]
  hgtEvent[[ GENE_ANNOTS ]] = hgtEvent[[ GENE_ANNOTS ]][ startPos:endPos ]
  return( hgtEvent )
}

# Find HGT events with redundant genes. If the genes for a given event
# are represented in a larger event, collapse them and remove the 
# smaller event from the output list
RemoveNestedEvents = function( hgtEventData )
{
  # Get the number of genes involved in each hgt event
  hgtGeneList = sapply( 1:length(hgtEventData), 
    function(i) list( hgtEventData[[i]][[ HGT_GENES ]] )
    )
  genesPerEvent = sapply( hgtGeneList, length )
  
  isNested = sapply( seq( length( hgtGeneList ) ), function(i)
  {
    toTest = which( genesPerEvent > genesPerEvent[i] )
    isNested = sapply( toTest, 
      function(j) !FALSE %in% ( hgtGeneList[[i]] %in% hgtGeneList[[j]] )
      )
    return( TRUE %in% isNested )
  })
  return( hgtEventData[ !isNested ] )
}

# ------------------------------------------------------------------------------