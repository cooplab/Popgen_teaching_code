#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)



simulate.pop<-function(N.vec=rep(5,30), const.RS=TRUE,  mutation= TRUE, mut.rate=  0.1, for.class= FALSE, initial.state="all.black",plot.freqs=FALSE,mult.pop=FALSE,pops=FALSE){
  cex.pch<- 0.9   #1.3
  #  c(rep(10,5),rep(3,2),rep(10,5),rep(3,2),rep(10,5))  #
  stopifnot(initial.state %in% c("all same","all different","two alleles","single mutation") )
  if(plot.freqs){ layout(c(1,2)); par(mar=c(1,2,0,1))}
  if(for.class){
    line.lwd<-1
    line.col<-"black"
    mut.line.lwd<-1
    mut.line.col<-"black"
    
  }else{
    line.lwd<-0.5
    line.col<-"grey"
    mut.line.lwd<-1
    mut.line.col<-"grey"
  }
  num.gens<- length(N.vec)-1	
  if(!mult.pop){
    ind.pop.par<-matrix(1,nrow=max(N.vec),ncol=num.gens+1)
    ind.pop<-matrix(1,nrow=max(N.vec),ncol=num.gens+1)
  }else{
    ind.pop.par<-pops[["ind.pop.par"]]
    ind.pop<-pops[["ind.pop"]]
  }
  
  num.gens<- length(N.vec)-1
  offset<-0.1
  plot(c(1,num.gens),c(0.5,max(N.vec))+c(-offset,offset),type="n",axes=FALSE,xlab="",ylab="")
  mtext(side=1,line=0,"Generations")
  text(1,0.5,"Past")
  text(num.gens-1,0.5,"Present")
  
  track.cols<- list()
  N <-N.vec[1]
  if(initial.state=="all same") my.cols<-rep("black",2*N)  #sample(rainbow(2*N))
  if(initial.state=="all different") my.cols<-sample(rainbow(2*N))
  if(initial.state=="two alleles")  my.cols<-  rep(c("blue","red"),N)
  if(initial.state=="single mutation")  my.cols<-  c("red",rep("blue",2*N-1))
  stopifnot((2*N)==length(my.cols))
  
  track.cols[[1]]<-my.cols
  points(rep(1,N),1:N+offset, pch=15,cex=cex.pch,col=my.cols[(1:N)*2])
  points(rep(1,N),1:N-offset, pch=15,cex=cex.pch,col=my.cols[(1:N)*2-1])
  
  for(i in 1:num.gens){
    
    N.new<-N.vec[i+1]
    N.old<-N.vec[i]
    points(rep(i,N.old),1:N.old+offset, pch=15,cex=cex.pch,col=my.cols[(1:N.old)*2])
    points(rep(i,N.old),1:N.old-offset, pch=15,cex=cex.pch,col=my.cols[(1:N.old)*2-1])
    new.cols<-rep("black",2*N.new)
    
    if(const.RS){ 
      repro.success<-rep(1/N.old,N.old)
    }else{
      repro.success<-sample(c(rep(0.5/(N.old),N.old-2),c(0.25,0.25)),replace=FALSE)
    }
    
    for(ind in 1:N.new){
      this.pop.par <- ind.pop.par[ind,i+1]
      available.pars <- (1:N.old)[which(ind.pop[1:N.old,i] == this.pop.par)]
      par<-sample(available.pars,2,replace=FALSE,prob=repro.success[which(ind.pop[1:N.old,i] == this.pop.par)])
      
      which.allele.1<-sample(c(-1,1),1)
      if(i != num.gens){ lines(c(i,i+1), c(par[1]+which.allele.1*offset,ind-offset),col=line.col,lwd=line.lwd)}
      new.cols[2*ind-1]<- my.cols[2*par[1] +ifelse(which.allele.1==1,0,-1)]
      
      which.allele.2<-sample(c(-1,1),1)
      if(i != num.gens){ lines(c(i,i+1), c(par[2]+which.allele.2*offset,ind+offset),col=line.col,lwd=line.lwd)}
      new.cols[2*ind]<- my.cols[2*par[2] +ifelse(which.allele.2==1,0,-1)]
      
      if(mutation){
        if(runif(1)<mut.rate){ 
          new.cols[2*ind-1]<- sample(rainbow(4*N),1)
          if(i != num.gens){ lines(c(i,i+1), c(par[1]+which.allele.1*offset,ind-offset),col=mut.line.col,lwd=mut.line.lwd)}
          
        }
        if(runif(1)<mut.rate){ 
          new.cols[2*ind]<- sample(rainbow(4*N),1)
          if(i != num.gens){ lines(c(i,i+1), c(par[2]+which.allele.2*offset,ind+offset),col=mut.line.col,lwd=mut.line.lwd)}
        } 
        
      }
    }	
    ##redraw points to cover lines		 
    points(rep(i,N.old),1:N.old+offset, pch=15,cex=cex.pch,col=my.cols[(1:N.old)*2])
    points(rep(i,N.old),1:N.old-offset, pch=15,cex=cex.pch,col=my.cols[(1:N.old)*2-1])
    my.cols<-new.cols
    track.cols[[i+1]]<-my.cols
    if(!const.RS) sapply(which(repro.success>1/N.old), function(ind){ draw.circle(x=i,y=ind,radius=0.2,nv=100,border=NULL,col=NA,lty=1,lwd=1)})
  }
  #	recover()
  if(plot.freqs){
    plot(c(1,num.gens),c(0,1),type="n",axes=FALSE,xlab="",ylab="")
    all.my.cols<-unique(unlist(track.cols))
    
    if(!mult.pop){ 
      my.col.freqs<-sapply(track.cols,function(my.gen){sapply(all.my.cols,function(my.col){sum(my.gen==my.col)})})
      sapply(all.my.cols,function(col.name){lines(my.col.freqs[col.name,]/(2*N.vec),col=col.name,lwd=2)});
    }else{
      
      for(pop in 1:max(ind.pop)){
        my.col.freqs<-sapply(1:num.gens, function(gen){
          #			recover()
          my.gen<-track.cols[[gen]]
          if(all(ind.pop.par[ind.pop[,gen]==pop,gen]==0)) return(rep(NA,length(all.my.cols)))  #if pop doesn't exist in this gen.
          
          these.inds<-which(ind.pop[,gen]==pop)
          my.gen<-c(my.gen[these.inds*2],my.gen[these.inds*2-1])
          sapply(all.my.cols,function(my.col){
            sum(my.gen==my.col)
          })})
        rownames(my.col.freqs)<-		all.my.cols
        sapply(all.my.cols[-length(all.my.cols)],function(col.name){lines(my.col.freqs[col.name,]/(2*5),col=col.name,lwd=2,lty=pop)});	
      }
    }
    
    axis(2)
  }
}

# Define UI for application that draws a histogram
#user interface
ui <- pageWithSidebar( 
  
  headerPanel = headerPanel("Genetic drift simulations"),
  
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
            selectInput("initialstate", "Initial state:",
                        choices = c("all same","all different","two alleles","single mutution")),
            sliderInput("mu",
                        "Mut. Rate mu:",
                        min = 0.0,
                        max = 1,
                        value = .01),
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
    #plot(c(0,1),c(0,1))
    #recover()
    #text(0.5,0.5,my.params$N)
    simulate.pop(N.vec=rep(my.params$N,20),initial.state=my.params$initialstate,
                               mut.rate=my.params$mu,plot.freqs=TRUE)   
   # simulate.pop(N.vec=rep(my.params$N,my.params$gens),initial.state=my.params$initial.state,
    #             mut.rate=my.params$mu,plot.freqs=TRUE)            
    
  })
  
  
}




# Run the application 
shinyApp(ui = ui, server = server)
