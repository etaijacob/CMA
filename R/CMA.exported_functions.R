#CMA.exported_functions.R

#' Correlated mutation analysis of amino acid and codon sequences by multiple methods
#'
#' \code{CMA_aa_codon} Calculates the Direct Coupling, OMES and MI measurments of all pairwise positions of
#' a protein sequence based on an input MSA file of codon and amino acid level sequences in fasta or sth format
#'
#' @param aa_alignment_file Input amino acid sequence alignment file name in a fasta or sth format
#' @param aa_alignment_file_type Amino acid alignment file format (fasta or sth)
#' @param codon_alignment_file Input codon sequence alignment file name in a fasta or sth format
#' @param codon_alignment_file_type Codon alignment file format (fasta or sth)
#' @param seqnameorindex Id or index of the sequence as given in the amino acid sequence alignment file to evalute the correlated mutations for
#' @param pseudocount_weight pseudo-count weight
#' @param theta theta
#' @return returns a data.frame of the measurements MI, OMES and DI at the amino and codon sequence level.
#' See below for a description of the measurements.
#'
#' @details
#' MI - **M**utual **I**nforamtion:
#' 1. Convert table MI data to a matrix.
#' 2. Execute the algorithm of MIp (below)
#' 3. Returns the same table with i & j and but the column, MIp
#'
#' Algorithm (based on http://bioinformatics.oxfordjournals.org/content/24/3/333.full):
#' MI(a, mean(x)) = 1/m(Sum(MI(a, x)))
#' where n is the number of columns in the alignment, m = n-1 for convenience,
#' and the summation is over x=1 to n, x != a.
#'
#' mean(MI) denotes the overall mean mutual information:
#' mean(MI) = 2/(m*n)(Sum(MI(x, y))),
#' where the indices run x=1 to m, y = x+1 to n.
#'
#' The average product correction:
#' APC(a, b) = (MI(a, mean(x))*MI(b, mean(x)))/mean(MI)
#' Is an approximation to the background MI shared by positions a and b.
#' MIp denotes the difference between total observed MI and the APC:
#' MIp(a, b) = MI(a, b) - APC(a, b).
#' In order to avoid negative values, values are shifted by adding the minimum
#'
#' @examples
#' #Performing correlated mutatios annalysis at the codon level
#' #for Q7K3Y9_DROME/661-715 given PF00014 full alignment:
#' res <- CMA_aa_codon(aa_alignment_file = system.file("extdata",
#'                                                     "PF00014_full.fasta.test",
#'                                                     package ="CMA"),
#'                     aa_alignment_file_type = "fasta",
#'                     codon_alignment_file = system.file("extdata",
#'                                                        "PF00014_full.fasta.test.tranalignout",
#'                                                        package ="CMA"),
#'                     codon_alignment_file_type = "fasta",
#'                     seqnameorindex = 1)
#' #Performing correlated mutatios annalysis for W5PT72_SHEEP/33-151 given
#' #PF00074 seed alignment
#' \dontrun{
#' res <- CMA_aa_codon(aa_alignment_file = system.file("extdata",
#'                                                     "PF00074_seed.sth",
#'                                                     package ="CMA"),
#'                     aa_alignment_file_type = "sth",
#'                     codon_alignment_file = system.file("extdata",
#'                                                        "PF00074_seed.sth.tranalignout",
#'                                                        package ="CMA"),
#'                     codon_alignment_file_type = "fasta",
#'                     seqnameorindex = 2)
#' }
#' @export
#'
CMA_aa_codon <- function(aa_alignment_file = system.file("extdata", "PF00014_full.fasta.test", package ="CMA"),
                         aa_alignment_file_type = "fasta",
                         codon_alignment_file = system.file("extdata", "PF00014_full.fasta.test.tranalignout", package ="CMA"),
                         codon_alignment_file_type = "fasta",
                         seqnameorindex = 1, pseudocount_weight = 0.5, theta = 0.2) {
  # 1. reads an sth file for a ref seq column indices
  # 2. dca aa file
  # 3. dca codon file
  if(aa_alignment_file_type == "fasta") {
    aaMSA <- getExcludedIndices.fasta(aa_alignment_file, seqname = seqnameorindex)
  } else if(aa_alignment_file_type == "sth") {
    aaMSA <- getExcludedIndices.sth(aa_alignment_file, seqnameorindex = seqnameorindex)
  }
  idxs <- aaMSA$idxs

  cat(sprintf("theta = %f, weight = %f.\n", theta, pseudocount_weight))
  cat("Doing AA..\n")
  dca.aa <- dca(inputfile = aa_alignment_file, seqid_or_excluded_indices = idxs,
                nuc = F, fileType = aa_alignment_file_type,
                pseudocount_weight = pseudocount_weight, theta = theta)
  cat("Doing CODON..\n")
  dca.nuc <- dca(inputfile = codon_alignment_file, seqid_or_excluded_indices = idxs,
                 nuc = T, fileType = codon_alignment_file_type,
                 pseudocount_weight = pseudocount_weight, theta = theta)
  res <- merge(dca.aa$results, dca.nuc$results, by = c("i","j"))
  names(res) <- c("i", "j", "MI.aa", "OMES.aa", "DI.aa", "MI.codon", "OMES.codon", "DI.codon")
  return(list(results=res, dca.aa=dca.aa, dca.codon=dca.nuc))
}


#' Correlated mutation analysis at the codon sequence level by multiple methods
#'
#' \code{CMA_codon} Calculates the Direct Coupling, OMES and MI measurments of all pairwise positions of
#' a protein sequence based on an input MSA file of codon level sequences in fasta or sth format
#'
#' @param aa_alignment_file Input amino acid sequence alignment file name in a fasta or sth format
#' @param aa_alignment_file_type Amino acid alignment file format (fasta or sth)
#' @param codon_alignment_file Input codon sequence alignment file name in a fasta or sth format
#' @param codon_alignment_file_type Codon alignment file format (fasta or sth)
#' @param seqnameorindex Id or index of the sequence as given in the amino acid sequence alignment file
#' to evalute the correlated mutations for
#' @return returns a data.frame of the measurements MI, OMES and DI at codon sequence level.
#' See below for a description of the measurements.
#'
#' @details
#' MI - **M**utual **I**nforamtion
#' calc MIp by:
#' 1. Convert table MI data to a matrix.
#' 2. Execute the algorithm of MIp (below)
#' 3. Returns the same table with i & j and but the column, MIp
#'
#' Algorithm (based on http://bioinformatics.oxfordjournals.org/content/24/3/333.full):
#' MI(a, mean(x)) = 1/m(Sum(MI(a, x)))
#' where n is the number of columns in the alignment, m = n-1 for convenience,
#' and the summation is over x=1 to n, x != a.
#'
#' mean(MI) denotes the overall mean mutual information:
#' mean(MI) = 2/(m*n)(Sum(MI(x, y))),
#' where the indices run x=1 to m, y = x+1 to n.
#'
#' The average product correction:
#' APC(a, b) = (MI(a, mean(x))*MI(b, mean(x)))/mean(MI)
#' Is an approximation to the background MI shared by positions a and b.
#' MIp denotes the difference between total observed MI and the APC:
#' MIp(a, b) = MI(a, b) - APC(a, b).
#' In order to avoid negative values, values are shifted by adding the minimum
#'
#' @examples
#' #Performing correlated mutatios annalysis at the codon level for Q7K3Y9_DROME/661-715 given
#' #PF00014 full alignment:
#' res.codon <- CMA_codon(aa_alignment_file = system.file("extdata",
#'                                                        "PF00014_full.fasta.test",
#'                                                        package ="CMA"),
#'                        aa_alignment_file_type = "fasta",
#'                        codon_alignment_file = system.file("extdata",
#'                                                           "PF00014_full.fasta.test.tranalignout",
#'                                                           package ="CMA"),
#'                        codon_alignment_file_type = "fasta",
#'                        seqnameorindex = "Q7K3Y9_DROME/661-715")
#' @export
CMA_codon <- function(aa_alignment_file = system.file("extdata", "PF00014_full.fasta.test", package ="CMA"),
                      aa_alignment_file_type = "fasta",
                      codon_alignment_file = system.file("extdata", "PF00014_full.fasta.test.tranalignout", package ="CMA"),
                      codon_alignment_file_type = "fasta",
                      seqnameorindex = "Q7K3Y9_DROME/661-715") {
  if(aa_alignment_file_type == "fasta") {
    aaMSA <- getExcludedIndices.fasta(aa_alignment_file, seqname = seqnameorindex)
  } else if(aa_alignment_file_type == "sth") {
    aaMSA <- getExcludedIndices.sth(aa_alignment_file, seqnameorindex = seqnameorindex)
  }
  idxs <- aaMSA$idxs
  cat("Doing AA..\n")
  dca.codon <- dca(codon_alignment_file, seqid_or_excluded_indices = idxs, nuc = T, fileType = codon_alignment_file_type)
  res <- dca.codon$results
  names(res) <- c("i", "j", "MI.codon", "OMES.codon", "DI.codon")
  return(list(results=res, dca.aa=dca.codon))
}

#' Correlated mutation analysis at the amino acid sequence level by multiple methods
#'
#' \code{CMA_codon} Calculates the Direct Coupling, OMES and MI measurments of all pairwise positions of
#' a protein sequence based on an input MSA file of codon level sequences in fasta or sth format
#'
#' @param aa_alignment_file Input amino acid sequence alignment file name in a fasta or sth format
#' @param aa_alignment_file_type Amino acid alignment file format (fasta or sth)
#' @param seqnameorindex Id or index of the sequence as given in the amino acid sequence alignment file
#' to evalute the correlated mutations for
#' @return returns a data.frame of the measurements MI, OMES and DI at amino acid sequence level.
#' See below for a description of the measurements.
#'
#' @details
#' MI - **M**utual **I**nforamtion:
#' 1. Convert table MI data to a matrix.
#' 2. Execute the algorithm of MIp (below)
#' 3. Returns the same table with i & j and but the column, MIp
#'
#' Algorithm (based on http://bioinformatics.oxfordjournals.org/content/24/3/333.full):
#' MI(a, mean(x)) = 1/m(Sum(MI(a, x)))
#' where n is the number of columns in the alignment, m = n-1 for convenience,
#' and the summation is over x=1 to n, x != a.
#'
#' mean(MI) denotes the overall mean mutual information:
#' mean(MI) = 2/(m*n)(Sum(MI(x, y))),
#' where the indices run x=1 to m, y = x+1 to n.
#'
#' The average product correction:
#' APC(a, b) = (MI(a, mean(x))*MI(b, mean(x)))/mean(MI)
#' Is an approximation to the background MI shared by positions a and b.
#' MIp denotes the difference between total observed MI and the APC:
#' MIp(a, b) = MI(a, b) - APC(a, b).
#' In order to avoid negative values, values are shifted by adding the minimum
#'
#' @examples
#' #Performing correlated mutatios annalysis at the
#' #amino acid level for Q7K3Y9_DROME/661-715 given PF00014 full alignment:
#' res.aa <- CMA_aa(aa_alignment_file = system.file("extdata",
#'                                                  "PF00014_full.fasta.test", package ="CMA"),
#'                  aa_alignment_file_type = "fasta",
#'                  seqnameorindex = "Q7K3Y9_DROME/661-715")
#' @export
CMA_aa <- function(aa_alignment_file = system.file("extdata", "PF00014_full.fasta.test", package ="CMA"),
                   aa_alignment_file_type = "fasta",
                   seqnameorindex = "Q7K3Y9_DROME/661-715") {
  if(aa_alignment_file_type == "fasta") {
    aaMSA <- getExcludedIndices.fasta(aa_alignment_file, seqname = seqnameorindex)
  } else if(aa_alignment_file_type == "sth") {
    aaMSA <- getExcludedIndices.sth(aa_alignment_file, seqnameorindex = seqnameorindex)
  }
  idxs <- aaMSA$idxs
  cat("Doing AA..\n")
  dca.aa <- dca(aa_alignment_file, seqid_or_excluded_indices = idxs, nuc = F, fileType = aa_alignment_file_type)
  res <- dca.aa$results
  names(res) <- c("i", "j", "MI.aa", "OMES.aa", "DI.aa")
  return(list(results=res, dca.aa=dca.aa))
}

