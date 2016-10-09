#' Query NCBI's dbSNP for information on a set of SNPs
#' 
#' @export
#' @param SNPs A vector of SNPs (rs numbers).
#' @param blocksize how many SNPs to look at each time
#' @param numcores cores for parallel
#' @param hownice how nice? options are: c("very", "somewhat", "not")
#' 

NCBI_snp_query_par <- function(SNPs, blocksize=200, numcores=36, hownice="very"){
  hownice <- match.arg(hownice, c("very", "somewhat", "not"))
  if(hownice=="very"){
    nice <- function(){Sys.sleep(runif(1, 1, 8))}
  } else if(hownice=="somewhat"){
    nice <- function(){Sys.sleep(runif(1, 1, 4))}
  } else if(hownice=="not") {
    nice <- function(){Sys.sleep(runif(1, 0.01, 0.5))}
  }
  
  GROUPS <- rep(1:(length(SNPs)/blocksize+1), blocksize)
  GROUPS <- sort(GROUPS[1:length(SNPs)])
  
  query_dbSNP <- function(group){
    nice()
    ret <- NCBI_snp_query2(SNPs[GROUPS==group])
    as.data.frame(ret$summary)
  }
  
  LDF <- mclapply(1:max(GROUPS),
	    query_dbSNP,
	    mc.cores=numcores)
  
  # mclapply returns try-error for failed threads
  
  is_not_DF <- sapply(LDF, class) != "data.frame"
  
  while(sum(is_not_DF) != 0) {
    print("some elements not data.frame, requerying")
    
    new_LDF <- mclapply( which(is_not_DF),
		  query_dbSNP,
		  mc.cores=numcores)
    
    for(s in 1:sum(is_not_DF)){
      LDF[[ which(is_not_DF)[s] ]] <- new_LDF[[s]]
    }
    
    is_not_DF <- sapply(LDF, class) != "data.frame"
  }
  
  do.call("rbind", LDF)
}
