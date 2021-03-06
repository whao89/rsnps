#' Get SNP (Single-Nucleotide Polymorphism) Data on the Web
#'
#' This package gives you access to data from OpenSNP (http://opensnp.org/) via their API
#' (https://opensnp.org/faq#api).
#'
#' @importFrom utils download.file read.table
#' @importFrom httr GET content stop_for_status
#' @importFrom plyr ldply llply laply compact
#' @importFrom stringr str_split str_replace_all str_trim
#' @importFrom XML xmlInternalTreeParse xmlToList 
#' @importFrom jsonlite fromJSON
#' @name rsnps-package
#' @aliases rsnps
#' @docType package
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
#' @author Kevin Ushey \email{kevinushey@@gmail.com}
NULL
