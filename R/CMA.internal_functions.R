#CMA.internal_functions.R

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
# In order to avoid negative values, values are shifted by adding the minimum
#Input: tbl: a data.frame with three columns; 1 - i, 2 - j, 3 - MI score
#            Other fields are ignored.
#Output:     same data.frame with column i & j and a column named "MIp"

calc_MIp <- function(tbl, eps = 0.00001) {

  tbl <- tbl[, 1:3]
  names(tbl) <- c("i", "j", "MI")
  tt <- rbind(tbl, data.frame(i=tbl$j, j=tbl$i, MI=tbl$MI))
  d <- reshape2::acast(tt, formula = i~j, add.missing=T)
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
  tbl <- tbl[, 1:3]
  names(tbl) <- c("i", "j", "MI")
  tt <- rbind(tbl, data.frame(i=tbl$j, j=tbl$i, MI=tbl$MI))
  d <- reshape2::acast(tt, formula = i~j, add.missing=T)
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
getExcludedIndices.sth <- function(sth_file, seqnameorindex) {
  sth <- read_stockholm_alignment(sth_file)
  if(is.na(sth[seqnameorindex, "seq"])) {
    cat(sprintf("ERROR - %s is not in %s.\n", seqnameorindex, sth_file))
    return(NA)
  }
  chars <- strsplit(sth[seqnameorindex,"seq"], "")[[1]]
  excluded.idxs <- which(grepl("[[:lower:]]+", chars) | chars == '.')

  return(list(sth=sth, idxs=excluded.idxs)) # last bug: which(!(grepl("[[:upper:]]+", chars) & chars != '.'))))
}


