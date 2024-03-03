
rm(list=ls())

package_list <- list("ggplot2", "dplyr", "tidyverse","devtools", "multcompView", "plotly", "streamgraph")

check_packages <- function(package) {
if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

lapply(package_list, check_packages)
lapply(package_list, library, character.only = TRUE)


# Install streamgraph from GitHub
#devtools::install_github("hrbrmstr/streamgraph")
##install.packages("streamgraph")
library(streamgraph)


source("functions.R")
