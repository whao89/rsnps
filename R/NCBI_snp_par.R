#' Query NCBI's dbSNP for information on a set of SNPs
#' 
#' @export
#' @param SNPs A vector of SNPs (rs numbers).
#' @param blocksize how many SNPs to look at each time
#' @param numcores cores for parallel
#' @param hownice how nice? options are: c("very", "somewhat", "not")
#' 

NCBI_snp_query_par <- function(SNPs, blocksize=25, numcores=28, hownice="very"){
  hownice <- match.arg(hownice, c("very", "somewhat", "not"))
  if(hownice=="very"){
    nice <- function(){Sys.sleep(runif(1, 3, 8))}
  } else if(hownice=="somewhat"){
    nice <- function(){Sys.sleep(runif(1, 0.25, 3))}
  } else if(hownice=="not") {
    nice <- function(){Sys.sleep(0.1)}
  }
  
  GROUPS <- rep(1:(length(SNPs)/blocksize+1), blocksize)
  GROUPS <- sort(GROUPS[1:length(SNPs)])
  LDF <- mclapply(1:max(GROUPS),
	    function(group){
	      nice()
	      ret <- NCBI_snp_query2(SNPs[GROUPS==group])
	      nice()
	      ret
	    },
	    mc.cores=numcores)
  do.call("rbind", LDF)
}
