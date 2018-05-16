library(shinydashboard)
library(beeswarm)
library(DT)
library(scater)
library(scde)
library(here)


scde_plot<-function(x, o.ifm, cd, o.prior){
  scde.test.gene.expression.difference(gene = x, models = o.ifm, counts = cd, prior = o.prior, return.details = F)
}
my_datatable<-function(x){datatable(x,
                                    extensions = 'Buttons',
                                    filter = 'top',
                                    selection = list(mode = "single"),
                                    rownames = TRUE,
                                    options = list(
                                      dom = "Bflrtip",
                                      lengthMenu = c(5,10,100,500,1000,2000,"All"),
                                      buttons =
                                        list("copy", list(
                                          extend = "collection",
                                          buttons = c("csv","excel","pdf","colvis"),
                                          text = "Download")
                                        ))
)}
cd <- readRDS('data/cd.rds')



sidebar <- dashboardSidebar(
  sidebarMenu(id="tabs",
              menuItem(
                "DGE",
                tabName = "dge",
                icon = icon("bar-chart-o")),
              conditionalPanel(
                "input.tabs == 'dge'",
                radioButtons("res1", "Comparisons", 
                             list("CS15/CS16 v CS20/FL8W",
                                  "CS15/CS16 v CB",
                                  "CB v CS20/FL8W"),
                             select = "CS15/CS16 v CS20/FL8W")
              )
  )
)


# body is dynamically rendered by server
body <- dashboardBody(tabItems(tabItem("dge",uiOutput("ui"))))
ui <- dashboardPage(
  dashboardHeader(title = "SC scde"),
  sidebar,
  body
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  rv <- reactiveValues(which_o.ifm=character(),
                       o.ifm=data.frame(),
                       which_ediff=character(),
                       ediff=data.frame(),
                       o.prior=character(),
                       results=character(),
                       gene="RUNX1")
  
observeEvent(input$res1,{ 
         if(input$res1=="CS15/CS16 v CS20/FL8W"){rv$which_o.ifm <- 'data/early_v_mid_o.ifm.rds'; 
                                                 rv$which_ediff <- 'data/early_vs_mid_ediff.rds'} 
    else if(input$res1=="CS15/CS16 v CB")       {rv$which_o.ifm <- 'data/late_v_early.o.ifm.rds';
                                                 rv$which_ediff <- 'data/early_vs_late_ediff.rds'}
    else if(input$res1=="CB v CS20/FL8W")       {rv$which_o.ifm <- 'data/late_vs_mid_o.ifm.rds';
                                                 rv$which_ediff <- 'data/late_vs_mid_ediff.rds'}
    rv$o.ifm <- readRDS(rv$which_o.ifm)
    rv$ediff <- readRDS(rv$which_ediff)
    rv$ediff$padj <- pnorm(-abs(rv$ediff$cZ))
    rv$ediff <- rv$ediff[, c("padj","ce")]
    rv$o.prior <- scde.expression.prior(models = rv$o.ifm, counts = cd, length.out = 400, show.plot = FALSE)


  })
  
  observeEvent(input$table_rows_selected,{rv$gene<-rownames(rv$ediff[input$table_rows_selected,])})
  output$table <- renderDataTable({my_datatable(rv$ediff)})
  
  output$diff_plot <- renderPlot({scde_plot(rv$gene, rv$o.ifm, cd, rv$o.prior )})
  
  output$ui <- renderUI({ fluidRow(box(title = input$res1, width = "80%", dataTableOutput("table")),
                                    box(title = rv$gene,plotOutput("diff_plot"), 
                                        downloadButton(outputId = "down", label = "Download the plot")))
    })
  
  
  # downloadHandler contains 2 arguments as functions, namely filename, content
  output$down <- downloadHandler(
    filename = paste(rv$gene,"png",sep="."),
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      png(file,width = 800, height = 600) # open the pdf device
      scde_plot(rv$gene, rv$o.ifm, cd, rv$o.prior ) # draw the plot
      dev.off()  # turn the device off
      
    } 
  )
   

}
# Run the application 
shinyApp(ui = ui, server = server)


