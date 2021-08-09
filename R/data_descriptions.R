#' A baseline community for examples
#'
#' The community contains seven nodes generally meant to represent oribatid mites feeding on microorganisms and detritus.
#' @format A community with a feeding matrix and properties dataframe:
#' \describe{
#'   \item{imat}{The feeding matrix.}
#'   \item{prop}{The properties data frame containing node names (ID), assimilation efficiency (a), production efficiency(p), C:N ratio (CN), biomass (B), death rate (d), proportion of death cycled back to a detrital pool (DetritusRecycling), Booleans stating whether the node is detritus, plant, and can immobilize nitrogen, and a list of mutual predators.}
#' }
#' @source \url{https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2435.13706}
"intro_comm"
