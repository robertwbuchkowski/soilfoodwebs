#' A function to compile the food web from simple data inputs
#'
#' @param  feeding A data frame listing the feeding relationships in the food web. Only needs to contain the nodes that have prey items along with their names. Three columns in this data frame must be: Predator, Prey, and Preference. The Predator column is the name of the predator, the Prey column is the name of the prey, and the Preference is the preference that the predator has for that prey item after correcting for abundance. Preference values are relative, so for each predator the preference value can be any number and all that matters is its proportion of the total of all preference values given to that predator. Preference is meaningless if a predator only has one prey item. A good default value is 1 for everything if you don't want to set preferences beyond prey abundance.
#' @param properties A data frame listing the properties or parameters in the food web. Must contain the names of all of the nodes in the food web. This data frame must contain: ID, d,a,p,B, CN, Detritusrecycling, isDetritus, isPlant, and canIMM.
#' @return A community that is compatible with the functions of soilfoodwebs. This is a list with the feeding matrix imat first and the properties prop second.
#' @examples
#'# Creating a simple three node community:
#'
#'# Create a data frame of feeding relationships:
#'feedinglist = data.frame(
#'Predator = c("Pred", "Pred"),
#'Prey = c("Prey1", "Prey2"),
#'Preference = c(1,1.2))
#'
#'# Create a data frame of properties for each species:
#'properties = data.frame(ID = c("Pred", "Prey1", "Prey2"), # Name
#'d = c(1,3,0.5), # Death rate
#'a = c(0.61,0.65,0.45), # Assimilation efficiency
#'p = c(0.5,0.4,0.3), # Production efficiency for carbon
#'B = c(0.1,8,5), # Biomass
#'CN = c(4.5,4.8,5), # Carbon to nitrogen ratio
#'DetritusRecycling = c(0,0,0), # proportion of detritus recycling
#'isDetritus = c(0,0,0), # Boolean: Is this pool detritus?
#'isPlant = c(0,0,0), # Boolean: Is this pool a plant?
#'canIMM = c(0,0,0)) # Boolean: Can the pool immobilize inorganic nitrogen?
#'
#'# Build the food web:
#'build_foodweb(feedinglist, properties)
#' @export

build_foodweb <- function(feeding,
                          properties){

  # Check that feeding has the right columns and values:
  if(!all(colnames(feeding) %in% c("Predator", "Prey", "Preference"))) stop("feeding needs to have the columns: Predator, Prey, and Preference")

  # Check that all of the feeding columns are the right type:
  if(!inherits(feeding$Preference, "numeric")) stop("feeding Preference needs to be numeric")
  if(!inherits(feeding$Predator, "character")) stop("feeding Preference needs to be a character")
  if(!inherits(feeding$Prey, "character")) stop("feeding Preference needs to be a character")

  # Check that there are na missing values:
  if(any(is.na(feeding$Preference))) stop("No feeding preference values can be NA")

  if(any(is.na((properties[,c("d","a","p","B", "CN", "DetritusRecycling", "isDetritus", "isPlant", "canIMM")])))) stop("No values in the properties data frame can be NA")

  # Check that all of the properties are present:
  if(!all(c("ID", "d","a","p","B", "CN", "DetritusRecycling", "isDetritus", "isPlant", "canIMM") %in% colnames(properties))) stop("The properties data frame must contain all of the following columns: ID, d,a,p,B, CN, DetritusRecycling, isDetritus, isPlant, and canIMM.")

  # Check that all of the feeding relationships are listed in the properties data frame:
  if(!(all(unique(c(feeding$Predator, feeding$Prey)) %in% properties$ID))) stop("Node names listed in 'feeding' are not present in the properties data frame as listed in the ID column. All nodes in the feeding matrix must have properties.")

  # Create the feeding matrix using the data from properties as a guide:

  feedingmatrix = matrix(0,
                         dimnames = list(properties$ID, properties$ID),
                         ncol = length(properties$ID), nrow = length(properties$ID), byrow = TRUE)

  #Run through the feeding list and add in the non-zero feeding links identified there.
  for(i in 1:dim(feeding)[1]){
    feedingmatrix[feeding$Predator[i],feeding$Prey[i]] = feeding$Preference[i]
  }

  # Build the community
  community <- list(imat = feedingmatrix,
                   prop = properties)

  # Check the community for accuracy
  community <- checkcomm(community)

  return(community)

}
