#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

library("RColorBrewer")

track_lineages<-function(N.vec, n.iter, num.tracked, col.allele=c("red","blue","purple"),return.tracked=FALSE){
    if(max(N.vec)<20){
        allele.col<-adjustcolor("black",0.7)
    }else{
        allele.col<- NA #adjustcolor("black",0.2)
    }
    if(num.tracked<9){ 
        col.allele<-brewer.pal(num.tracked,"Dark2")
    }else{
        col.allele<-rainbow(num.tracked)
    }    
    offset<-0.2
    num.gens<-length(N.vec)
    for(iter in 1:n.iter){
        N.max<-max(N.vec)
        N<-N.vec[num.gens]
        N.prev<-N.vec[num.gens-1]
        plot(c(1,num.gens),c(1,N.max),type="n",axes=FALSE,xlab="",ylab="")
        mtext(side=1,line=1,"Generations")
        
        track.this.allele<-vector("list", 2*N)
        track.this.allele.time<-list()
        track.this.allele[sample(1:(2*N),num.tracked)]<-1:num.tracked
        
        track.this.allele.next.gen<-vector("list", 2*N.prev)
        
        for(i in num.gens:2){
            if(return.tracked) track.this.allele.time[[i]]<-track.this.allele
            N<-N.vec[i]
            N.prev<-N.vec[i-1]
            track.this.allele.next.gen<-vector("list", 2*N.prev)
            for(ind in 1:N){
                
                par<-sample(1:N.prev,2,replace=FALSE)
                which.allele<-sample(c(-1,1),1)
                lines(c(i,i-1), c(ind-offset,par[1]+which.allele*offset),col="light grey",lwd=0.5)
                if(!is.null(track.this.allele[[2*ind-1]])){
                    this.one<-2*par[1] +ifelse(which.allele==1,0,-1); 
                    track.this.allele.next.gen[[this.one]]  <- c(track.this.allele.next.gen[[this.one]],track.this.allele[[2*ind-1]])
                }
                
                which.allele<-sample(c(-1,1),1)
                lines(c(i,i-1), c(ind+offset,par[2]+which.allele*offset),col="light grey",lwd=0.5)
                if(!is.null(track.this.allele[[2*ind]])){ 
                    this.one<-2*par[2] +ifelse(which.allele==1,0,-1); 
                    track.this.allele.next.gen[[ this.one]]  <- c(track.this.allele.next.gen[[this.one]],track.this.allele[[2*ind]])
                }
                #		recover()
            }
            for(this.allele in 1:num.tracked){ 
                daughter<-which(sapply(track.this.allele,function(allele){any(allele==this.allele)}))
                parent<-which(sapply(track.this.allele.next.gen,function(allele){any(allele==this.allele)}))
                lines(c(i,i-1), c(ceiling(daughter/2)+offset* ifelse(daughter %% 2,-1,1) ,ceiling(parent/2) + offset*ifelse(parent %% 2,-1,1) ),col=col.allele[this.allele],lwd=2)
            }
            
            points(rep(i,N),1:N+offset, pch=19,cex=1,col=allele.col)
            points(rep(i,N),1:N-offset, pch=19,cex=1,col=allele.col)
            track.this.allele<-track.this.allele.next.gen
        }
        
        
    }
    #recover()
    if(return.tracked) track.this.allele.time
}


# Define UI for application that draws a histogram
#user interface
ui <- pageWithSidebar( 
    
    headerPanel = headerPanel("Coalescent simulations"),
    
    sidebarPanel(
        # Sidebar with a slider input for number of bins 
        sliderInput("gens",
                    "Num. gens:",
                    min = 2,
                    max = 100,
                    value = 30),
        sliderInput("N",
                    "Pop size N:",
                    min = 1,
                    max = 100,
                    value = 10),
        sliderInput("n",
                    "Sample size n:",
                    min = 1,
                    max = 100,
                    value = 3),
 #       selectInput("initialstate", "Initial state:",
 #                   choices = c("all same","all different","two alleles","single mutution")),
        actionButton("goButton", "GO"),
        #plotOutput(outputId = "freq"),
        #downloadButton(outputId = 'drift_out.csv', label = 'Download')
    ), 
    mainPanel =  mainPanel(
        plotOutput(outputId = 'freq')
    )
)

#back end code and response to user input
server <- function(input, output){
    
    #each time user hits "go"
    rand <- eventReactive(input$goButton, {
        #parameters
        
        #package data for plotting
        return(input)
    })
    
    
    output$freq <- renderPlot({
        my.params<- rand()
    #    plot(c(0,1),c(0,1))
        #recover()
    #    text(0.5,0.5,my.params$N)
        track_lineages(N.vec=rep(my.params$N,my.params$gens), n.iter=1, num.tracked=my.params$n)
        
        #simulate.pop(N.vec=rep(my.params$N,20),initial.state=my.params$initialstate,
         #            mut.rate=my.params$mu,plot.freqs=TRUE)   
        # simulate.pop(N.vec=rep(my.params$N,my.params$gens),initial.state=my.params$initial.state,
        #             mut.rate=my.params$mu,plot.freqs=TRUE)            
        
    })
    
    
}



# Run the application 
shinyApp(ui = ui, server = server)
