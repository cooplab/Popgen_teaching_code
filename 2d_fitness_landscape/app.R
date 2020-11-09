#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(mvtnorm)


x<-seq(0,20,length=200)
y<-seq(0,20,length=200)



norm.surf<-function(my.means,my.sigma){
    norm.2D<-outer(x,y,FUN=function(x,y){dmvnorm(x=cbind(x,y),mean=my.means,sigma = my.sigma*diag(length(my.means)))})
    print(max(norm.2D))
    norm.2D
}


plot_landscape<-function(){
    my.surf.1<-norm.surf(c(15,7),my.sigma=5)
    my.surf.2<-norm.surf(c(7,12),my.sigma=5)
    my.surf<- my.surf.1 + my.surf.2*1.2
    image(x,y,my.surf,xlim=c(2,18),ylim=c(2,18),xlab="",ylab="",axes=FALSE,cex.lab=1.8,asp=1) 
    mtext("Phenotype 1",side=1,cex=1.8)
    mtext("Phenotype 2",side=2,cex=1.8)
    return(my.surf)
}
    

plot_ascent<-function(my_xs,my_ys,my.surf,arrow.col="black",genetic_corr){
    Vx<-100
    Vy<-100
    Vxy<-genetic_corr*sqrt(Vy * Vx)
    points(my_xs,my_ys)
    if(!is.null(my_xs)){ #text(my_x[1],my_y[1],"blah")
        for(i in 1:length(my_xs)){    
            
            my_x=my_xs[i]; my_y=my_ys[i]
            for(gen in 1:50){
                x_int<-which(my_x >= x[-length(x)] & my_x < x[-1])
                y_int<-which(my_y >= y[-length(y)] & my_y < y[-1])
                
                beta_x <- (my.surf[x_int+1,y_int] - my.surf[x_int,y_int])/(x[x_int+1]-x[x_int])
                beta_y <- (my.surf[x_int,y_int+1] - my.surf[x_int,y_int])/(y[y_int+1]-y[y_int])
                
                R_x <- beta_x*Vx + beta_y*Vxy
                R_y <- beta_y*Vy + beta_x*Vxy
                
                new.dist<-sqrt(R_x^2+R_y^2)
                x_new<- my_x + R_x
                y_new<- my_y + R_y
                #arrows(x0=my_x,x1=my_x+1,y0=my_y,y1=my_y+1,length=0.12*1,lwd=1.2,col=arrow.col) 
                arrows(x0=my_x,x1=x_new,y0=my_y,y1=y_new,length=0.12*new.dist,lwd=1.2,col=arrow.col)  #if(abs(R_x)>0.03 |abs(R_y)>0.03)
                my_x<-x_new
                my_y<-y_new
                
            }
        }
    }
}

add.genetic.corr.points<-function(genetic_corr){
    
    lines(c(2,6),c(4,4))
    lines(c(4,4),c(4-2,4+2))
    
    sd.matrix<-0.4*diag(2)
    sd.matrix[1,2]<- genetic_corr * 0.4
    sd.matrix[2,1]<- sd.matrix[1,2]
    my.norms<-rmvnorm(500,mean=c(4,4),sigma=sd.matrix);
    points(my.norms[,1],my.norms[,2],col=adjustcolor("black",0.2),pch=19)
    
}

ui <- fluidPage(
    fluidRow(
               plotOutput("plot1", click = "plot_click")
        ),
    fluidRow(
    sidebarPanel(
        sliderInput("genetic_corr",
                    "Genetic Correlation:",
                    min = -0.99,
                    max = 0.99,
                    value = 0)
    ),
    actionButton("updateplot", "Update Plot:")
    )
    
)

# ui <- pageWithSidebar( 
#     
#     headerPanel = headerPanel("2D selection landscape"),
#     
#     sidebarPanel(
#         sliderInput("genetic_corr",
#                     "Genetic Correlation:",
#                     min = -1,
#                     max = 1,
#                     value = 0)
#     ),
#     mainPanel(
#         plotOutput("plot1", click = "plot_click"), 
#     )
# )

#ui <- basicPage(
#    plotOutput("plot1", click = "plot_click"),
    #verbatimTextOutput("info"),
#    actionButton("updateplot", "Update Plot:")
#)

server <- function(input, output) {
    val <- reactiveValues(clickx = NULL, clicky = NULL)
    
    observe({
        input$plot_click
        isolate({
            val$clickx =  c(val$clickx, input$plot_click$x)
            val$clicky = c(val$clicky, input$plot_click$y)     
        })
    }) #adding clicks to list
    
    rand <- eventReactive(input$goButton, {
        #parameters
        val<-numeric()
        #package data for plotting
        return(input)
    })
    
    
    output$plot1 <- renderPlot({
        my.surf<-plot_landscape() #plot(mtcars$wt, mtcars$mpg)
        add.genetic.corr.points(genetic_corr=input$genetic_corr)
        plot_ascent(val$clickx, val$clicky,my.surf=my.surf,genetic_corr=input$genetic_corr)
        
        
 #       input$updateplot
  #      isolate({
   #         val<-numeric()
    #    })
        #isolate({
           #points(val$clickx, val$clicky)
        #})
    })
    
    output$info <- renderText({
        paste0("x = ", val$clickx, ", y = ",val$clicky, "\n")
    })
    
}

shinyApp(ui, server)

# 
# ui <- basicPage(
#     plotOutput("plot1", click = "plot_click"),
#     verbatimTextOutput("info")
# )
# 
# server <- function(input, output) {
#     output$plot1 <- renderPlot({
# #        plot_landscape()
# #        if(!is.null(input$plot_click$x)) plot.ascent(input$plot_click$x,input$plot_click$y)
#     })
#     
#     #output$addTraj<- plot.ascent(input$plot_click$x,input$plot_click$y)
#         
#     output$info <- renderText({
#         paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
#     })
# }
# 
# shinyApp(ui, server)
