#!/usr/bin/Rscript
## Author: Mahek Pannu
## mkpannu@bu.edu
## BU BF591
## Final Project 

library("tidyverse")


load_counts <- function(filename) {
  counts <- read.table(filename) %>%
    rownames_to_column("gene") %>%
    as_tibble()
  return(counts)
}