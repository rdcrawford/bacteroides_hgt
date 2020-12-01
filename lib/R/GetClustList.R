# ------------------------------------------------------------------------------
# GetClustList
# 2018/07/20
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Higherarcical clustering is applied to a matrix and the matrix is then 
# divided into k-clusters
# ------------------------------------------------------------------------------

GetClustList = function(
  inputMatrix, # Matrix with row names to be list elemnets in the output list
  kVal,        # Number of clusters to divide the matrix into
  borderCols   # Optional. Vecotr with colors for the output plot
  )
{
  # ---- Parse the input arguments ---------------------------------------------
  
  if ( missing(borderCols) )
  {
    # Assign a random color to each cluster
    borderCols = rainbow( kVal )[ sample.int( kVal, kVal) ]
    
  } else if ( length(borderCols) != kVal ) {
    
    stop( "The number of colors to be plotted must be the same as the kVal" )
  }
  
  # ---- Function Definitions --------------------------------------------------

  hClustFunc = function(x) hclust( x, method = "complete" )
  distFunc   = function(x) as.dist( ( 1 - cor( t(x) ) ) / 2 )

  # ---- Divide the input matrix into clusters ---------------------------------
  
  # Create a distance matrix and apply apply the clustering algorithm 
  distMat  = distFunc( inputMatrix )
  clustMat = hClustFunc( distMat )
  
  # Divide into the requested number of clusters
  clusts = cutree( clustMat, k = kVal )
  
  # Plot the dendrogram and draw a box around each cluster
  plot( clustMat )
  rect.hclust( clustMat, k = kVal, border = borderCols )
  
  # Create a list of the representitives within each cluster
  return( tapply( names( clusts ), clusts, function(x) return(x) ) )
}

# ------------------------------------------------------------------------------

