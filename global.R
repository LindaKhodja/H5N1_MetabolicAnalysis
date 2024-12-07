#global.R

library(shiny)
library(shinydashboard)
library(visNetwork)
library(xml2)
library(RColorBrewer)
library(DT)

#fonctions utilitaires
source("utils.R")

# Palette de couleurs pour les compartiments
compartment_colors <- c(
  "cytosol" = "darkred",
  "extracellular region" = "darkgreen",
  "endoplasmic reticulum membrane" = "darkorange",
  "endoplasmic reticulum lumen" = "magenta",
  "Golgi membrane" = "gold",
  "endocytic vesicle membrane" = "yellowgreen",
  "nucleoplasm" = "pink",
  "endosome lumen" = "cyan",
  "plasma membrane" = "skyblue",
  "unknown" = "blue",
  "nuclear envelope" = "purple",
  "mitochondrial inner membrane" = "brown",
  "mitochondrial intermembrane space" = "gray"
)