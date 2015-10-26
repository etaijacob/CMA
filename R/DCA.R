#DCA.R


#' Direct Coupling Analysis (DCA) for amino acid and codon sequences

#'   @param  inputfile file containing the FASTA alignment
#'   @param pseudocount_weight - relative weight of pseudo count
#'   @param theta threshold for sequence id in reweighting
#'   @param seqid Three options:
#'                 1) A single numeric value indicates the sequence index in the MSA to use as a reference
#'                   (1 is the default - An MSA column is either all dot+lower case or dash+upper case,
#'                    by very definition of the output of HMMer.
#'                 2) A string indicates the name or id of the sequence in the msa to use as a reference
#'                 3) A vector of intergers indicates the columns to exclude from analysis
#'                    (no use of reference sequence in this case).
#'   @param nuc TRUE for codon based analysis and FALSE (default) for amino acids based analysis.
#'
#'  @return a list containing data produced in this function and the dca results.
#'          The results is composed by N(N-1)/2
#'          (N = length of the sequences) rows and 4 columns:
#'          residue i (column 1), residue j (column 2),
#'          MI(i,j) (Mutual Information between i and j), and
#'          DI(i,j) (Direct Information between i and j).
#'          Note: all insert columns are removed from the alignment.
#'   @seealso \code{\link[MAPDB]}
#'   @seealso \url{https://github.com/etaijacob/AA2CODON} if you want to generate the input files
#'
#' @examples
#' cma.res <- dca(codon_file_name, seqid = 1, nuc = T, fileType = "fasta")
#'
#'  SOME RELEVANT VARIABLES:
#'   N        number of residues in each sequence (no insert)
#'   M        number of sequences in the alignment
#'   Meff     effective number of sequences after reweighting
#'   q        equal to 21 (20 aminoacids + 1 gap)
#'   align    M x N matrix containing the alignmnent
#'   Pij_true N x N x q x q matrix containing the reweigthed frequency
#'            counts.
#'   Pij      N x N x q x q matrix containing the reweighted frequency
#'            counts with pseudo counts.
#'   C        N(q-1) x N(q-1) matrix containing the covariance matrix.
#'
#########################################################################
#' @export
dca <- function(inputfile, pseudocount_weight = 0.5, theta = 0.2,
                seqid = 1, nuc = F, fileType="fasta") {

  cat("return_alignment.\n")
  print(system.time(msa <- return_alignment(inputfile, seqid = seqid, nuc = nuc, fileType = fileType)))
  cat("Compute_True_Frequencies.\n")
  print(system.time(tfs <-  Compute_True_Frequencies(msa$Z, msa$M, msa$N, msa$q, theta)))
  cat(sprintf("### N = %d M = %d Meff = %.2f q = %d\n", msa$N, msa$M, tfs$Meff, msa$q))
  cat("with_pc\n")
  print(system.time(pc <- with_pc(tfs$Pij_true, tfs$Pi_true, pseudocount_weight, msa$N, msa$q)))
  cat("Compute_C\n")
  print(system.time(C <- Compute_C(pc$Pij, pc$Pi, msa$N, msa$q)))
  cat("InvC.\n")
  #print(system.time(invC <- solve(C)))
  #inverse via the Choleski decomposition
  print(system.time(invC <- chol2inv(chol(C))))
  cat("Compute_Results.\n")
  print(system.time(results <- Compute_Results(pc$Pij, pc$Pi, tfs$Pij_true, tfs$Pi_true, invC, msa$N, msa$q, outputfile)))

  return(list(msa = msa,
              tfs = tfs,
              pc = pc,
              C = C,
              invC = invC,
              results = results))
}


#all characters should be an uppercase and AAs or Nucs only
#mind that this function treats aa differently than nuc file:
#in nuc file it excludes all gaps '---' and '...' and 'XXX' where in aa file only '.'.
#Best thing is to give as input the sequence indices to exclude in seqid parameter.
#' @export
return_alignment <- function(inputfile="Y://PFAM2/data/sth/seed/PF00001_v27_seed.sth.aligndprots",
                             seqid = 1, nuc = F, fileType = c("fasta", "sth"), doUnique=T) {
  # reads alignment from inputfile, removes inserts and converts into numbers
  cat("return_alignment..\n")
  #require(seqinr)
  attributes(seqid) <- list(type="Reference sequence index")
  if(fileType == "fasta") {
    msa <- read.alignment(inputfile, format = "fasta", forceToLower = F)
    cat(sprintf("Original MSA length is: %d.\n", length(msa$seq)))
  } else if(fileType == "sth") {
    msa <- read.stockholm.alignment(inputfile)
    cat(sprintf("Original MSA length is: %d.\n", length(msa$seq)))
  }

  columnidxs <- NULL
  if(length(seqid) > 1 & is.numeric(seqid)) {
    cat("seqid is a numeric vector\n")
    columnidxs <- seqid
    attributes(seqid) <- list(type="Excluded columns")

  } else if(is.numeric(seqid)) {
    if(seqid > length(msa$nam)) {
      cat(sprintf("%d is greater than %d. Taking the first as ref.\n", seqid, length(msa$nam)))
      seqid <- 1
    }
  } else {
    print(seqid)
    myid <- grep(seqid, msa$nam)[1]

    if(length(myid) == 0) {
      cat(sprintf("%s does not exists in MSA. Taking the first seq as ref.\n", seqid))
      seqid <- 1
    } else {
      cat(sprintf("Ref id found: %d out of %d possibilities\n", myid, length(grep(seqid, msa$nam))))
      seqid <- myid
      if(doUnique) {
        if(fileType == "fasta") {
          msa$nam  <- c(msa$nam[seqid], msa$nam)
          msa$seq  <- c(msa$seq[seqid], msa$seq)
          msa$nb  <- c(msa$nb[seqid], msa$nb)
          msa$com  <- c(msa$com[seqid], msa$com)
        } else if(fileType == "sth") {
          msa$name  <- c(msa$name[seqid], msa$name)
          msa$seq  <- c(msa$seq[seqid], msa$seq)
        }
        seqid <- 1
      }
    }
  }
  if(doUnique) {
    if(fileType == "sth") {
      msa <- msa[!duplicated(msa$seq),]
    } else if(fileType == "fasta") {
      inidxs <- which(!duplicated(msa$seq))
      msa$nb  <- msa$nb[inidxs]
      msa$nam <- msa$nam[inidxs]
      msa$seq <- msa$seq[inidxs]
      msa$com <- msa$com[inidxs]
    }
  }
  cat(sprintf("After redundancy (unique) MSA length is: %d.\n", length(msa$seq)))
  cat(sprintf("Ref sequence is %s after removal of duplicates.\n", paste(seqid, collapse = ",")))
  if(nuc == FALSE) {
    tmp <- do.call(rbind, lapply(strsplit((msa$seq), ""), function(x) factor(x, levels = aas)))
    print(dim(tmp))
    #return(tmp)
    if(is.null(columnidxs)) {
      if(fileType == "sth") { #excluding only '.' and lower case aas from analysis as defined by HMMER
        excludeidxs <- which(tmp[seqid, ] == 1 | is.na(tmp[seqid, ]))
      } else if(fileType == "fasta") { #exclude gaps of ref sequence from analysis
        excludeidxs <- which(tmp[seqid, ] == 22)
        print(tmp[seqid,])
      }
      if(length(excludeidxs) > 0)
        tmp <- tmp[, -excludeidxs] #exclude gaps/dots from analysis
    } else {
      if(max(columnidxs) > dim(tmp)[2]) {
        cat(sprintf("ERROR - there are only %d columns where input column idx max val is %s.\n", dim(tmp)[2], max(columnidxs)))
        return(NA)
      }
      tmp <- tmp[, -columnidxs] #exclude gaps (in hammer means '.' and lower case letters) from analysis
    }
    tmp[ tmp > 21 | is.na(tmp) ] <- 1 # 'X' and '-' turn to be 1

  } else {
    tmp <- do.call(rbind, strsplit((msa$seq), ""))
    im <- matrix(1:dim(tmp)[2], ncol=dim(tmp)[2]/3)
    tmp <- do.call(rbind, lapply(1:dim(tmp)[1],
                                 function(z) sapply(lapply(1:dim(im)[2],
                                                           function(x) tmp[z, ][im[,x]]),
                                                    function(x) factor(paste(x, collapse=""), levels=codons))))
    if(is.null(columnidxs)) {
      tmp <- tmp[, -which(tmp[seqid, ] == 1 | is.na(tmp[seqid, ]))] #exclude gaps from analysis
    } else {
      if(max(columnidxs) > dim(tmp)[2]) {
        cat(sprintf("ERROR - there are only %d columns where input column idx max val is %s.\n", dim(tmp)[2], max(columnidxs)))
        return(NA)
      }
      tmp <- tmp[, -columnidxs] #exclude gaps from analysis

    }
    print(class(tmp))
    tmp[ tmp > 65 | is.na(tmp) ] <- 1 # 'xxx' and '...' turn to be 1
  }
  rownames(tmp) <- msa$nam
  align_full <- tmp
  M <- dim(align_full)[1]

  N <- dim(align_full)[2]
  q <- max(align_full)
  return(list(N=N, M=M, q=q, Z=align_full, seqid = seqid, raw.msa = msa))
}


# Two different implementations:
#   1. One for RW
#   TODO - 2. One without - further optimizations can be done (using table)
Compute_True_Frequencies <- function(align,M,N,q,theta, useSampling=F) {
  # computes reweighted frequency counts
  # the number of nonredundant sequences is measured as the effective sequence number Meff after reweighting (see Methods).
  # The comparison to results
  # without reweighting and to reweighting at 80%
  cat("Compute_True_Frequencies..\n")
  W <- rep(1, M)

  if( theta > 0.0 ) {
    if(useSampling) {
      W <- 1/(1 + colSums(sampledHammingdist(align, Ns = 240) < theta))
    } else {
      W <- 1/(1 + colSums(hammingdist(align) < theta))
    }
  }
  Meff <- sum(W)

  cat(sprintf("Allocating: %d, %d withb Meff = %f\n", N, q, Meff))

  #   Pi_true <- array(0, c(N, q))

  #    for(j in 1:M)
  #      for(i in 1:N)
  #        Pi_true[i, align[j,i]] <- Pi_true[i,align[j,i]] + W[j]

  Pi_true <- AccPi_true(align, W, N, q, M)
  Pi_true <- Pi_true/Meff

  cat("Going for main loop..\n")
  # Pij_true <- array(0, c(N, N, q, q))

  #   for(l in 1:M)
  #     for(i in 1:(N-1))
  #       for(j in (i+1):N) {
  #         Pij_true[i, j, align[l,i],align[l,j]] <- Pij_true[i, j, align[l, i], align[l, j]] + W[l]
  #         Pij_true[j, i, align[l,j],align[l,i]] <- Pij_true[i, j, align[l, i], align[l, j]]
  #       }

  Pij_true = AccPij_true(align, W, N, q, M)
  Pij_true <- Pij_true/Meff

  scra <- diag(q)
  for(i in 1:N)
    for(alpha in 1:q)
      for(beta in 1:q)
        Pij_true[i, i, alpha, beta] <- Pi_true[i, alpha] * scra[alpha, beta]

  return(list(Pij_true=Pij_true, Pi_true= Pi_true, Meff = Meff, W = W))
}



#Sampled Hamming distance, the percentage of coordinates that differ in a sample of Ns comparisons.
#Input X: matrix with n rows vectors
#Output m: each column j represents a Ns sampled distances from seq j
sampledHammingdist <- function(X, Ns) {
  cat("sampledHammingdist..\n")
  n <- nrow(X)
  if(n < Ns) {
    print("is less")
    m <- matrix(nrow=n, ncol=n)
    for(i in seq_len(n))
      for(j in seq(i, n))
        m[j, i] <- m[i, j] <- sum(X[i,] != X[j,])
    return(m/ncol(X))
  } else {
    print("is bigger")
    getmydist <- function(i) {
      sapply(sample(seq(1, n)[-i], Ns), function(j) sum(X[i,] != X[j,]))
    }
    m <- sapply(seq_len(n), getmydist)/ncol(X)
    return(m)
  }
}


with_pc <- function(Pij_true, Pi_true, pseudocount_weight, N, q) {
  # adds pseudocount
  cat("with_pc..\n")
  Pij <- (1 - pseudocount_weight) * Pij_true + pseudocount_weight / q / q * array(1, c(N, N, q, q))
  Pi <- (1 - pseudocount_weight) * Pi_true + pseudocount_weight / q * array(1, c(N, q))

  scra <- diag(q)
  for(i in 1:N)
    Pij[i, i, , ] <- (1 - pseudocount_weight) * Pij_true[i, i, ,] + pseudocount_weight / q * scra

  return(list(Pij = Pij, Pi = Pi))
}


mapkey <- function(i, alpha, q) {
  A <- (q-1) * (i-1) + alpha
  return(A)
}




ReturnW <- function(C, i, j, q) {
  # extracts coupling matrix for columns i and j
  #cat("ReturnW..\n")
  W <- array(1, c(q,q))
  W[1:(q-1), 1:(q-1)] <- exp( -C[mapkey(i, 1:(q-1), q), mapkey(j, 1:(q-1), q)] )
  #print(W)
  return(W)
}

bp_link <- function(i, j, W, P1, q) {
  # computes direct information
  #cat("bp_link..\n")
  mus <- compute_mu(i, j, W, P1, q)
  DI <- compute_di(i, j, W, mus$mu1, mus$mu2, P1)
  return(DI)
}


compute_mu <- function(i, j, W, P1, q) {
  #cat("compute_mu..\n")
  epsilon <- 1e-4
  diff <- 1.0
  mu1 <- array(1,q)/q
  mu2 <- array(1,q)/q
  pi <- P1[i, ]
  pj <- P1[j, ]

  while ( diff > epsilon ) {
    scra1 <- as.numeric(mu2 %*% t(W))
    scra2 <- as.numeric(mu1 %*% W)

    new1 <- as.vector(pi / scra1)
    new1 <- as.vector(new1 / sum(new1))

    new2 <- as.vector(pj / scra2)
    new2 <- as.vector(new2 / sum(new2))

    diff <- max( abs( new1 - mu1 ), abs( new2 - mu2 ) )
    #print(diff)
    mu1 <- new1
    mu2 <- new2
  }
  #cat(sprintf("%d, %d: %f %f\n", i, j, mu1, mu2))
  return(list(mu1 = mu1, mu2 = mu2))
}


compute_di <- function(i, j, W, mu1,mu2, Pia) {
  # computes direct information
  #cat("compute_di..\n")
  tiny <- 1.0e-100
  Pdir <- W * ((as.matrix(mu1)) %*% t(as.matrix(mu2))) #direction of vector is the opposite than matlab

  #Pdir <- W * as.vector(t(mu1) %*% mu2)
  #print(Pdir)
  Pdir <- Pdir / sum(Pdir)
  #print(Pdir)
  #print(dim(Pdir))
  Pfac <- as.matrix(Pia[i,]) %*% t(as.matrix(Pia[j, ])) #direction of vector is the opposite than matlab
  #print(Pfac)
  #cat(sprintf("%d, %d: %f %f %f\n", i, j, Pdir, Pfac, (Pdir+tiny)/(Pfac+tiny)))
  DI <- sum(diag( t(Pdir) %*% log( (Pdir+tiny)/(Pfac+tiny) ) ))
  #print(log(Pdir+tiny)/(Pfac+tiny))
  return(DI)
}


Compute_Results <- function(Pij, Pi, Pij_true, Pi_true, invC, N, q, fnameout = NULL) {
  # computes and prints the mutual and direct informations
  cat("Compute_Results..\n")
  df <- NULL

  nns <- combn(N, 2)
  df <- do.call(rbind,
                lapply(1:dim(nns)[2],
                       function(x) {
                         MI <- calculate_mi_and_omes(nns[1, x] - 1, nns[2, x] - 1, Pij_true, Pi_true, q, N)
                         W_mf <- ReturnW(invC, nns[1, x], nns[2, x], q)
                         DI_mf_pc <- bp_link(nns[1, x], nns[2, x], W_mf, Pi, q)
                         return(data.frame(i = nns[1, x], j = nns[2, x], MI = MI[1], OMES = MI[2], DI = DI_mf_pc))
                       }))

  cat(sprintf("Results dim: (%d, %d)", dim(df)[1], dim(df)[2]))
  if(!is.null(fnameout)) {
    cat(sprintf("Writing output to file: %s.\n", fnameout))
    write.table(df, file = fnameout, sep = "\t", row.names = F, quote = F)
  }
  return(df)
}


#reads a stockholm file format
# only lines which do not start with a # are considered
# returns a dataframe with name and seq
#' @export
read.stockholm.alignment <- function(file="Y://PFAM2/data/sth/seed/PF00013_v27_seed.sth") {
  mm <- readLines(file)
  mm <- do.call(rbind, strsplit(mm[-grep("^#", mm)], "\\s+", perl = T))
  mm <- data.frame(name=mm[,1], seq=mm[,2], row.names = mm[,1], stringsAsFactors = F)
  return(mm)
}

codon_msa_to_binary_matrix <- function(nmat) {
  tt <- lapply(1:dim(nmat)[1], function(x) codon_seq_to_binary_vec(nmat[x,]))
  tt <- do.call(rbind, tt)
  rownames(tt) <- rownames(nmat)
  return(tt)
}

codon_seq_to_binary_vec <- function(sequence) {
  unlist(lapply(sequence, codon_letter_to_binary_vec))
}

codon_letter_to_binary_vec <- function(codon) {
  (1:65==codon)+0
}


aa_msa_to_binary_matrix <- function(pmat) {
  tt <- lapply(1:dim(pmat)[1], function(x) aa_seq_to_binary_vec(pmat[x,]))
  tt <- do.call(rbind, tt)
  rownames(tt) <- rownames(pmat)
  return(tt)
}

aa_seq_to_binary_vec <- function(sequence) {
  unlist(lapply(sequence, aa_letter_to_binary_vec))
}

aa_letter_to_binary_vec <- function(aa) {
  (1:21==aa)+0
}


aas <- toupper(c(".","A", "C", "D", "E", "F", "G", "H", "I",
                 "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W",
                 "Y", "-", "X"))

aas.v2 <- toupper(c("-","A", "C", "D", "E", "F", "G", "H", "I",
                    "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W",
                    "Y", ".", "X"))

codons <-toupper( c("---", "ttt", "ttc", "tta", "ttg", "tct", "tcc", "tca", "tcg", "tat", "tac",
                    "taa", "tag", "tgt", "tgc", "tga", "tgg", "ctt", "ctc", "cta",
                    "ctg", "cct", "ccc", "cca", "ccg", "cat", "cac", "caa", "cag",
                    "cgt", "cgc", "cga", "cgg", "att", "atc", "ata", "atg", "act",
                    "acc", "aca", "acg", "aat", "aac", "aaa", "aag", "agt", "agc",
                    "aga", "agg", "gtt", "gtc", "gta", "gtg", "gct", "gcc", "gca",
                    "gcg", "gat", "gac", "gaa", "gag", "ggt", "ggc", "gga", "ggg", "...", "xxx"))

