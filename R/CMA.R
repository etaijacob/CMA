#CMA.R

#' CMA: An R package for correlated mutation analysis both at codon and amino acid levels.
#'
#' @details
#' CMA provides several methods for calculating correlated mutations analysis.
#' Methods that are currently provided in this package are: OMES, MI, DCA.
#'
#' If you use this package please cite:
#'
#' Etai Jacob, Ron Unger and Amnon Horovitz (2015). Codon-level information improves predictions of inter-residue contacts
#' in proteins by correlated mutation analysis.
#' eLife 2015;10.7554/eLife.08932 URL \url{http://dx.doi.org/10.7554/eLife.08932}.
#'
#' It has several functions for correlated mutation analysis:
#' \itemize{
#' \item Correlated mutation analysis at the codon and amino acid levels
#' \item Read fasta and stockholm format files
#' }
#'
#' To learn more about CMA, start with the vignettes:
#' \code{browseVignettes(package = "CMA")}
#'
#' @seealso \code{\link[MAPDB]{PDBs2Pfam}} if you want to calculate the accuracy of your predictions based on solved structures.
#'
#' @seealso \url{https://github.com/etaijacob/AA2CODON} if you want to generate the input files
#'
#' Some notes about MSA characters and their use:
#'
#' HMM models often contain dots (".") and lower-case letters, which are preserved in the aligned-alignment view.
#' Alignment columns with dots or lower-case letters are not present in the HMM profile (see the Pfam FAQ).
#' Thus, an HMM profile can contain less columns than its seed alignment.
#' From PFAM FAQ \url{http://pfam.xfam.org/help#tabview=tab3}
#' The '-' and '.' characters both represent gap characters.
#' However they do tell you some extra information about how the HMM has generated the alignment.
#' The '-' symbols are where the alignment of the sequence has used a delete state in the HMM to jump past a match state.
#' This means that the sequence is missing a column that the HMM was expecting to be there.
#' The '.' character is used to pad gaps where one sequence in the alignment has sequence from the HMMs insert state.
#' See the alignment below where both characters are used. The HMM states emitting each column are shown.
#' Note that residues emitted from the Insert (I) state are in lower case.
#'
#' FLPA_METMA/1-193     ---MPEIRQLSEGIFEVTKD.KKQLSTLNLDPGKVVYGEKLISVEGDE
#'
#' FBRL_XENLA/86-317    RKVIVEPHR-HEGIFICRGK.EDALVTKNLVPGESVYGEKRISVEDGE
#'
#' FBRL_MOUSE/90-321    KNVMVEPHR-HEGVFICRGK.EDALFTKNLVPGESVYGEKRVSISEGD
#'
#' O75259/81-312        KNVMVEPHR-HEGVFICRGK.EDALVTKNLVPGESVYGEKRVSISEGD
#'
#' FBRL_SCHPO/71-303    AKVIIEPHR-HAGVFIARGK.EDLLVTRNLVPGESVYNEKRISVDSPD
#'
#' O15647/71-301        GKVIVVPHR-FPGVYLLKGK.SDILVTKNLVPGESVYGEKRYEVMTED
#'
#' FBRL_TETTH/64-294 KTIIVK-HR-LEGVFICKGQ.LEALVTKNFFPGESVYNEKRMSVEENG
#'
#' FBRL_LEIMA/57-291 AKVIVEPHMLHPGVFISKAK.TDSLCTLNMVPGISVYGEKRIELGATQ
#'
#' Q9ZSE3/38-276 SAVVVEPHKVHAGIFVSRGKsEDSLATLNLVPGVSVYGEKRVQTETTD
#'
#' HMM STATES MMMMMMMMMMMMMMMMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMM
#'
#' In this program:
#' Lower-case letters and (".") characters are filtered out, upper-case letters and "-" characters are preserved.
#' Other characters are not allowed. Only the 20 basic IUPAC protein single-letters are allowed.
#'
#'
#' This implementation of DCA is based mainly on Andrea and Martin's MATLAB source code (See below).
#'
#'
#' Copyright for this R implementation:
#' 2015 - Etai Jacob \email{etai.jacob@@gmail.com}
#'
#' Copyright for the matlab implementation which the DCA R source code is based on:
#' 2011/12 - Andrea Pagnani  \email{andrea.pagnani@@gmail.com} and Martin Weigt \email{martin.weigt@@upmc.fr}
#'
#'
#'
#' Further description of the DCA algorithm can be found at:
#'
#'
#'   F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander,
#'   R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
#'   analysis of residue co-evolution captures native contacts across
#'   many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
#'
#'
#' @docType package
#' @name CMA
#' @useDynLib CMA
#' @importFrom Rcpp sourceCpp
#' @import seqinr
#'
NULL
