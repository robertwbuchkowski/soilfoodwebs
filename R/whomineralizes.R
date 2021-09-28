#' Direct and indirect contributions to mineralizations
#'
#' @param usin The community in which we want to calculate mineralization rates.
#' @return A table of node effects on mineralization rates.
#' @details
#' The results are labeled as follows with direct contributions calculated from the full food web and indirect contributions calculated from the food web without that node. Indirect contributions do not include the direct contribution (i.e., it is subtracted).
#' * DirectC: The direct contribution to carbon mineralization.
#' * DirectN: The direct contribution to nitrogen mineralization.
#' * IndirectC: The direct contribution to carbon mineralization.
#' * IndirectN: The direct contribution to nitrogen mineralization.
#'
#' @examples
#' # Basic example for the introductory community:
#' whomineralizes(intro_comm)
#' @export
whomineralizes <- function(usin){
  usin = checkcomm(usin)
  Nnodes = dim(usin$imat)[1]
  Nnames = usin$prop$ID

  res1 <- comana(usin)

  output = data.frame(
    ID = unname(usin$prop$ID),
    DirectC = res1$Cmin/sum(res1$Cmin),
    DirectN = rowSums(res1$Nmin)/sum(res1$Nmin),
    IndirectC = NA,
    IndirectN = NA)
  rownames(output) = NULL
  for(rmnode in Nnames){
    usinmod = removenodes(usin, rmnode)
    res2 = comana(usinmod)
    output[output$ID == rmnode, "IndirectC"] =
      (sum(res1$Cmin) - sum(res2$Cmin) - output[output$ID == rmnode, "DirectC"])/sum(res1$Cmin)
    output[output$ID == rmnode, "IndirectN"] =
      (sum(res1$Nmin) - sum(res2$Nmin) - output[output$ID == rmnode, "DirectN"])/sum(res1$Nmin)
  }
  return(output)
}

