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
