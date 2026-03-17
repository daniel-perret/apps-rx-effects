#### libraries and setup code

library(tidyverse)
library(rFIA)
library(ggplot2)
library(terra)
library(sf)
library(ggmap)
#library(rgdal)
library(lme4)
library(performance)
library(ggeffects)
library(plotrix)
library(tidyterra)
library(RSQLite)

select <- dplyr::select

# setting my preferred ggplot2 theme
theme_set(theme_bw())
theme_update(text = element_text(size=16, color = "black"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             panel.border=element_rect(linewidth=1.5))

options(scipen = 9999)
#gcinfo(TRUE)
