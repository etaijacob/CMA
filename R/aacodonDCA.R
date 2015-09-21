#aacodonDCA.R


#' @export
codonDCA.fasta <- function(aa_fasta_file = "data/PF00014_full.fasta.test",
                           codon_file = "data/PF00014_full.fasta.test.tranalignout",
                           seqname = "Q7K3Y9_DROME/661-715") {
  fastaMSA <- getExcludedIndices.fasta(aa_fasta_file, seqname = seqname)
  idxs <- fastaMSA$idxs
  cat("Doing AA..\n")
  dca.codon <- dca(codon_file, seqid = idxs, nuc = T, fileType = "fasta")
  res <- dca.codon$results
  names(res) <- c("i", "j", "MI.codon", "OMES.codon", "DI.codon")
  return(list(results=res, dca.aa=dca.codon))
}

#' @export
aaDCA.fasta <- function(fasta_file = "data/PF00014_full.fasta.test",
                        seqname = "Q7K3Y9_DROME/661-715") {
  fastaMSA <- getExcludedIndices.fasta(fasta_file, seqname = seqname)
  idxs <- fastaMSA$idxs
  cat("Doing AA..\n")
  dca.aa <- dca(fasta_file, seqid = idxs, nuc = F, fileType = "fasta")
  res <- dca.aa$results
  names(res) <- c("i", "j", "MI.aa", "OMES.aa", "DI.aa")
  return(list(results=res, dca.aa=dca.aa))
}

#' @export
codonDCA.sth <- function(sth_file = "data/PF00074_seed.sth", codon_file= "data/PF00074_seed.sth.tranalignout",
                         codon_file_type = "fasta", seqname = 1) {
  sthMSA <- getExcludedIndices.sth(sth_file, seqname = seqname)
  idxs <- sthMSA$idxs
  cat("Doing AA..\n")
  dca.codon <- dca(codon_file, seqid = idxs, nuc = T, fileType = codon_file_type)
  res <- dca.codon$results
  names(res) <- c("i", "j", "MI.codon", "OMES.codon", "DI.codon")
  return(list(results=res, dca.codon=dca.codon))
}

#' @export
aaDCA.sth <- function(sth_file = "data/PF00074_seed.sth", seqname = 1) {
  sthMSA <- getExcludedIndices.sth(sth_file, seqname = seqname)
  idxs <- sthMSA$idxs
  cat("Doing AA..\n")
  dca.aa <- dca(sth_file, seqid = idxs, nuc = F, fileType = "sth")
  res <- dca.aa$results
  names(res) <- c("i", "j", "MI.aa", "OMES.aa", "DI.aa")
  return(list(results=res, dca.aa=dca.aa))
}


#seqnameorindex can be an integer (e.g. 1) or a seq name (e.g. "FMR1_DROME/286-359")
#In Hmmer MSA it is equivalent (see DCA comment header)
#' @export
aacodonDCA <- function(sth_file = "data/PF00074_seed.sth",
                       aa_file = "data/PF00074_seed.sth.aligndprots",
                       codon_file = "data/PF00074_seed.sth.tranalignout",
                       seqnameorindex=1, pseudocount_weight = 0.5, theta = 0.2, idxs = NULL) {
  # 1. reads an sth file for a ref seq column indices
  # 2. dca aa file
  # 3. dca codon file

  if(is.null(idxs)) {
    sthMSA <- getExcludedIndices.sth(sth_file, seqnameorindex)
    idxs <- sthMSA$idxs
  }
  cat(sprintf("theta = %f, weight = %f.\n", theta, pseudocount_weight))
  cat("Doing AA..\n")
  dca.aa <- dca(inputfile = aa_file, seqid = idxs, nuc = F,
                pseudocount_weight = pseudocount_weight, theta = theta)
  cat("Doing CODON..\n")
  dca.nuc <- dca(inputfile = codon_file, seqid = idxs, nuc = T,
                 pseudocount_weight = pseudocount_weight, theta = theta)
  res <- merge(dca.aa$results, dca.nuc$results, by = c("i","j"))
  names(res) <- c("i", "j", "MI.aa", "OMES.aa", "DI.aa", "MI.codon", "OMES.codon", "DI.codon")
  return(list(results=res, dca.aa=dca.aa, dca.codon=dca.nuc))
}

#calc MIp by:
#1. Convert table MI data to a matrix.
#2. Execute the algorithm of MIp (below)
#3. Returns the same table with i & j and but the column, MIp
#
#Algorithm (based on http://bioinformatics.oxfordjournals.org/content/24/3/333.full):
#MI(a, mean(x)) = 1/m(Sum(MI(a, x)))
#where n is the number of columns in the alignment, m = n-1 for convenience,
#and the summation is over x=1 to n, x != a.
#
#mean(MI) denotes the overall mean mutual information:
#mean(MI) = 2/(m*n)(Sum(MI(x, y))),
#where the indices run x=1 to m, y = x+1 to n.
#
#The average product correction:
#APC(a, b) = (MI(a, mean(x))*MI(b, mean(x)))/mean(MI)
#Is an approximation to the background MI shared by positions a and b.
#MIp denotes the difference between total observed MI and the APC:
#MIp(a, b) = MI(a, b) - APC(a, b).
#
# In order to avoid negative values values are shifted by adding the minimum
#Input: tbl: a data.frame with three columns; 1 - i, 2 - j, 3 - MI score
#            Other fields are ignored.
#Output:     same data.frame with column i & j and a column named "MIp"

calc_MIp <- function(tbl, eps = 0.00001) {
  require(reshape2)
  tbl <- tbl[, 1:3]
  names(tbl) <- c("i", "j", "MI")
  tt <- rbind(tbl, data.frame(i=tbl$j, j=tbl$i, MI=tbl$MI))
  d <- acast(tt, formula = i~j, add.missing=T)
  n = dim(d)[1]
  meanMI = sapply(1:n, function(x) mean(d[x,], na.rm = T))
  overallMI = mean(d[(upper.tri(d))])

  getMyMIp <- function(x) {
    APC <- (meanMI[tbl$i[x]]*meanMI[tbl$j[x]])/overallMI
    MIp <- tbl$MI[x] - APC
    return(MIp)
  }

  res <- data.frame(i = tbl$i, j = tbl$j, MIp = sapply(1:dim(tbl)[1], getMyMIp))
  #fix values < 0:
  res$MIp <- res$MIp + abs(min(res$MIp)) + eps

  return(res)

}

calc_APC <- function(tbl) {
  require(reshape2)
  tbl <- tbl[, 1:3]
  names(tbl) <- c("i", "j", "MI")
  tt <- rbind(tbl, data.frame(i=tbl$j, j=tbl$i, MI=tbl$MI))
  d <- acast(tt, formula = i~j, add.missing=T)
  n = dim(d)[1]
  meanMI = sapply(1:n, function(x) mean(d[x,], na.rm = T))
  overallMI = mean(d[(upper.tri(d))])

  getMyAPC <- function(x) {
    APC <- (meanMI[tbl$i[x]]*meanMI[tbl$j[x]])/overallMI
    return(APC)
  }

  res <- data.frame(i = tbl$i, j = tbl$j, APC = sapply(1:dim(tbl)[1], getMyAPC))
  return(res)

}

#For fast file only
#Exclude '-' from reference seq
#' @export
getExcludedIndices.fasta <- function(fasta_file,
                                     seqname) {
  msa <- seqinr::read.alignment(fasta_file, format = "fasta", forceToLower = F)
  fasta <- data.frame(name=msa$nam, seq=msa$seq, stringsAsFactors = F)

  idx <- grep(seqname, fasta$name)
  cat(sprintf("Found seq name %s as index %d.\n", seqname, idx))
  if(length(idx) == 0) {
    cat(sprintf("ERROR - %s is not in %s.\n", seqname, fasta_file))
    return(NA)
  }
  #Taking only the first occurance:
  chars <- strsplit(fasta[idx,"seq"], "")[[1]]

  excluded.idxs <- which(chars == '-')

  return(list(fasta=fasta, idxs=excluded.idxs)) # last bug: which(!(grepl("[[:upper:]]+", chars) & chars != '.'))))
}

#for sth file only:
#Exclude '.' and lower case letters
#' @export
getExcludedIndices.sth <- function(sth_file, seqnameorindex) {
  sth <- read.stockholm.alignment(sth_file)
  if(is.na(sth[seqnameorindex, "seq"])) {
    cat(sprintf("ERROR - %s is not in %s.\n", seqnameorindex, sth_file))
    return(NA)
  }
  chars <- strsplit(sth[seqnameorindex,"seq"], "")[[1]]
  excluded.idxs <- which(grepl("[[:lower:]]+", chars) | chars == '.')

  return(list(sth=sth, idxs=excluded.idxs)) # last bug: which(!(grepl("[[:upper:]]+", chars) & chars != '.'))))
}

