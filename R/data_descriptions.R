#' A baseline community for examples
#'
#' The community contains seven nodes generally meant to represent oribatid mites feeding on microorganisms and detritus.
#' @format A community with a feeding matrix and properties dataframe:
#' \describe{
#'   \item{imat}{The feeding matrix. Rows eat columns.}
#'   \item{prop}{The properties data frame containing node names (ID), assimilation efficiency (a), production efficiency(p), C:N ratio (CN), biomass (B), death rate (d), proportion of death cycled back to a detrital pool (DetritusRecycling), Booleans stating whether the node is detritus, plant, and can immobilize nitrogen, and a list of mutual predators.}
#' }
#' @source \url{https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2435.13706}
"intro_comm"

#' The soil food web published for CPER
#'
#' The community contains 17 nodes in bacterial, fungal, and root energy channels.
#'
#' @format A community with a feeding matrix and properties dataframe:
#' \describe{
#'   \item{imat}{The feeding matrix. Rows eat columns.}
#'   \item{prop}{The properties data frame containing node names (ID), assimilation efficiency (a), production efficiency(p), C:N ratio (CN), biomass (B), death rate (d), proportion of death cycled back to a detrital pool (DetritusRecycling), Booleans stating whether the node is detritus, plant, and can immobilize nitrogen, and a list of mutual predators. Biomass is in kilograms of carbon per hectare and turnover/death rate is in years.}
#' }
#' @source \url{https://doi.org/10.1007/BF00260580}
"CPER"

#' The soil food webs published for grazed and ungrazed plots in the Shortgrass Steppe long-term research station.
#'
#' The community contains 21 nodes in bacterial, fungal, and root energy channels.
#'
#' @format A list of six communities each with a feeding matrix and properties dataframe. NOTE: You must select one of the communities from the list to use the package functions. For example: comana(Andres2016$GA) carries out calculations for the GA plot.
#' \describe{
#'   \item{GA}{The grazed plot A.}
#'   \item{UGA}{The ungrazed plot A.}
#'   \item{GB}{The grazed plot B.}
#'   \item{UGB}{The ungrazed plot B.}
#'   \item{GC}{The grazed plot C.}
#'   \item{UGC}{The ungrazed plot C.}
#'   \item{imat}{Within the community: The feeding matrix. Rows eat columns.}
#'   \item{prop}{Within the community: The properties data frame containing node names (ID), assimilation efficiency (a), production efficiency(p), C:N ratio (CN), biomass (B), death rate (d), proportion of death cycled back to a detrital pool (DetritusRecycling), Booleans stating whether the node is detritus, plant, and can immobilize nitrogen, and a list of mutual predators. Biomass is in kilograms of carbon per hectare and turnover/death rate is in years.}
#' }
#'
#' @source \url{https://doi.org/10.1016/j.soilbio.2016.02.014}
"Andres2016"

#' The soil food web published for CPER
#'
#' The community contains 41 nodes in bacterial, fungal, and root energy channels. One node has zero biomass and is removed by the check_comm() function.
#'
#' @format A community with a feeding matrix and properties dataframe:
#' \describe{
#'   \item{imat}{The feeding matrix. Rows eat columns.}
#'   \item{prop}{The properties data frame containing node names (ID), assimilation efficiency (a), production efficiency(p), C:N ratio (CN), biomass (B), death rate (d), proportion of death cycled back to a detrital pool (DetritusRecycling), Booleans stating whether the node is detritus, plant, and can immobilize nitrogen, and a list of mutual predators. Biomass is in milligrams of carbon per square meter and turnover/death rate is in years.}
#' }
#' @source \url{http://link.springer.com/10.1007/s00300-017-2201-5}
"Koltz2018"

#' The soil food webs published along a chronosequence in the Netherlands.
#'
#' The community contains 21 nodes in bacterial, fungal, and root energy channels.
#'
#' @format A list of four communities each with a feeding matrix and properties dataframe. NOTE: You must select one of the communities from the list to use the package functions. For example: comana(Holtkamp2011$Young) carries out calculations for the Young field.
#' \describe{
#'   \item{Young}{The Young field.}
#'   \item{Mid}{The middle age field.}
#'   \item{Old}{The old field.}
#'   \item{Health}{The health.}
#'   \item{imat}{Within the community: The feeding matrix. Rows eat columns.}
#'   \item{prop}{Within the community: The properties data frame containing node names (ID), assimilation efficiency (a), production efficiency(p), C:N ratio (CN), biomass (B), death rate (d), proportion of death cycled back to a detrital pool (DetritusRecycling), Booleans stating whether the node is detritus, plant, and can immobilize nitrogen, and a list of mutual predators. Biomass is in kilograms of carbon per hectare to 10-cm depth and turnover/death rate is in years.}
#' }
#'
#' @source \url{https://doi.org/10.1016/j.soilbio.2010.10.004}
"Holtkamp2011"
