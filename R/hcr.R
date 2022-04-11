#' @title Apply management stratey (MS)
#
#' @description Apply management strategy (MS)
#'
#' @param F HERE:
#' @param MStype HERE:
#' @param SSB HERE:
#' @param Bpa HERE:
#'
#' @export
apply.ms <- function(F, MStype = "noF", SSB = NA, Bpa = NA){
    switch(as.character(MStype),
           ## Fish at FMSY
           "noF" = {
               return(0.0)
           },
           ## Fish at FMSY
           "FMSY" = {
               return(F)
           },
           ## Fish at FMSY with linear reduction if SSB < Bpa
           "ICES" = {
               return(ifelse(SSB > Bpa, F, F * SSB/Bpa))
           },
           stop("MStype not known!")
           )
}
