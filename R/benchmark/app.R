#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(TMB)
library(tmbstan)
library(cmdstanr)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Benchmarking some stats packages"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("params",
                        "Number of params:",
                        min = 1,
                        max = 50,
                        value = 30),
            actionButton("button","Run the model")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    randomVals <- eventReactive(input$button,{
      source('run_benchmark_shiny.R')
      run_models(input$params)
      setWinProgressBar(pb, i/(100)*100, label=info)
      })
                 output$distPlot <- renderPlot({
      #withProgress(message = 'Running models', value = 0, {
        # generate bins based on input$bins from ui.R
        
        #n <- 1
        #incProgress(1/n, detail = paste("Doing part", i))})
                   hist(randomVals())
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
