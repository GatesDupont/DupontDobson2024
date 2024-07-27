library(tidyverse)
library(jsonlite)
library(vroom)


# ---- Functions ----

# Get taxnomy CSV file from eBird API
get_taxonomy <- function(){
  
  url <- "https://api.ebird.org/v2/ref/taxonomy/ebird"
  
  data <- data.table::fread(url) %>%
    select(`Common Name` = COMMON_NAME, 
           species_code = SPECIES_CODE)
  
}

# Get taxnomy CSV file from eBird API
get_sci_names <- function(){
  
  url <- "https://api.ebird.org/v2/ref/taxonomy/ebird"
  
  data <- data.table::fread(url) %>%
    select(`Common Name` = COMMON_NAME, 
           sci_name = SCIENTIFIC_NAME)
  
}
