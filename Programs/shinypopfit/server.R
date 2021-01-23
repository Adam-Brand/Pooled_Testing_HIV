#==============================================================================
# FILENAME: server.R
# PROJECT: 	Pooled testing in HIV
# PURPOSE: this is the server file for the shinypopfit shiny application
#          
#          
# AUTHOR: Adam Brand


# INPUT datasets: none

# OUPUT: none - this program runs a shiny app locally


# R VERSION: 3.6.1
#==============================================================================
#Notes: to run, click the "Run App" button on the top right of an R Studio session with this file loaded





# =============================================================================


library(shiny)
source("shinysource.R")

shinyServer(function(input, output){
  
  output$plot2 <- renderPlot({
    set.seed(input$seed)
    
    plot.pop(n=input$popsize, prev.over.cutoff=input$prev, cutoff=input$cutoff, 
             shape1=input$shape1, shape2=input$shape2, pf.prev=input$pf, b0=input$b0, 
             b1=input$b1, b2=input$b2, b3=input$b3, sd=input$sd)
  })
})
