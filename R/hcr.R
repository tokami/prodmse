#' @title Apply harvest control rule (HCR)
#
#' @description Apply harvest control rule (HCR)
#'
#' @param HCR HCR type
#' @param biomass Curent biomass (either SSB or ESB) should correspond to biomass type in refs!
#' @param refs Reference points, should contain Fmsy, Bpa, Blim
#'
#' @export
apply.hcr <- function(HCR = "noF", biomass = NULL, biomass.blim = NULL, refs = NULL){
    switch(as.character(HCR),
           ## Fish at FMSY
           "noF" = {
               return(0.0)
           },
           ## Fish at FMSY
           "FMSY" = {
               return(refs$Fmsy)
           },
           ## Fish at FMSY with linear reduction if SSB < Btrigger
           "ICES" = {
               return(ifelse(biomass > refs$Btrigger, refs$Fmsy, refs$Fmsy * biomass/refs$Btrigger))
           },
           "ICES2" = {
               if(is.null(biomass.blim)){
                   return(ifelse(biomass > refs$Btrigger, refs$Fmsy, refs$Fmsy * biomass/refs$Btrigger))
               }else{
                   return(ifelse(biomass.blim < refs$Blim, 0,
                          ifelse(biomass > refs$Btrigger, refs$Fmsy, refs$Fmsy * biomass/refs$Btrigger)))
               }
           },
           "ICES3" = {
               return(ifelse(biomass > refs$Btrigger, refs$Fmsy, refs$Fmsy * biomass/refs$Btrigger))
           },
           stop("HCR not implemented!")
           )
}
