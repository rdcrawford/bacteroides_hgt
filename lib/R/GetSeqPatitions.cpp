#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// MultiSeqAlgn
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
std::list< std::vector< std::string > >  GetSeqPatitions(
  std::vector< std::string > seqs, 
  std::vector< int >         genePartitions
  )
{
  // Initialize the list, reserving the number of matricies to output
  std::list< std::vector< std::string > > partitionList;

  // Initialize the start of the gene partition
  unsigned int gStart;
  unsigned int numSeqs = seqs.size();
  unsigned int len;

  // For each set of patitions
  for ( unsigned  int i = 0; i < genePartitions.size(); i++ )
  {
    if ( i == 0 ) gStart = 0;
    else gStart = genePartitions[ i - 1 ];
    len = genePartitions[ i ] - gStart;

    // Create an alignment with only the sequences in the partitions
    std::vector< std::string > subSeqs( numSeqs );
    for ( unsigned int j = 0; j < numSeqs; j++ )
      subSeqs[ j ] = seqs[ j ].substr( gStart, len );

    // Set the row and column names of the matirx and Return
    partitionList.push_back( subSeqs );
  }

  return partitionList;
}

// -----------------------------------------------------------------------------