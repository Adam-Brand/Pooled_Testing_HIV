#==============================================================================
# FILENAME: ui.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: ui.R file for the shiny program shinypopfit
#          
#          
# AUTHOR: Adam Brand


# INPUT datasets: none

# OUTPUT: none, this is a user interface for the a shiny app


# R VERSION: 3.6.1
#==============================================================================
#Notes: this needs to be in the same folder as the serverR file for this app, but does not need to be run





# =============================================================================



# This code also needs to be in a separate file, but called ui.R
library(shiny)
source("shinysource.R")

shinyUI(fluidPage(
  titlePanel("Parametric Population Fit"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("seed", label = "Seed", value=12, min=1),
      numericInput("popsize", label = "Pop. Size", value=1000, min=1000),
      numericInput("b0", label = "B0", value=1, min=-50),
      numericInput("b1", label = "B1", value=.05, min=-50),
      numericInput("b2", label = "B2", value=1, min=-50),
      numericInput("b3", label = "B3", value=.05, min=-50),
      numericInput("sd", label = "Standard Deviation", value=.5, min=0.01),
      sliderInput("cutoff", label = "Failure Cutoff", min=200, max=2500, value=1000, step=50),
      sliderInput("prev", label = "Prevalence Over Cutoff", min=.01, max=1, value=0.1),
      sliderInput("pf", label = "Prior Failure Prevalence", min=.01, max=1, value=0.25),
      numericInput("shape1", label = "Beta Shape 1", min=0.1, value=5),
      numericInput("shape2", label = "Beta Shape 2", min=0.1, value=0.5)
    ),
    mainPanel(
      h1("Population vs Parametric Fit: Combined curves"),
      plotOutput("plot2")
    )
  )
))
