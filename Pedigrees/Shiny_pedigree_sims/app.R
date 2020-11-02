#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#   Shiny server and panel setup borrowed from Silas Tittes https://github.com/silastittes/shiny_popgen

library(shiny)
library("kinship2")
library("RColorBrewer")


plot.ped.allele<-function(my.ped,which.allele,allele.cols,my.cex=3){
    num<-length(allele.cols)
    if(which.allele==1) my.ped$affected=data.frame(allele1=rep(1,num),allele2=rep(0,num))
    if(which.allele==2) my.ped$affected=data.frame(allele1=rep(0,num),allele2=rep(1,num))
    if(which.allele==0) my.ped$affected=data.frame(allele1=rep(0,num),allele2=rep(0,num))
    
    my.ped.obj <- pedigree(
        id=my.ped$id, 
        dadid=my.ped$dadid, 
        momid=my.ped$mumid,
        sex=my.ped$sex, 
        famid=my.ped$ped,
        affected=as.matrix(my.ped$affected)
    )  
    if(which.allele %in% c(0,2)) par(new=TRUE)
    plot(my.ped.obj['100'],cex=my.cex,id=rep(NA,num),col=allele.cols,density=c(-1,-1))  #c("kid1","kid2","dad","mum")
    
}


generate.alleles<-function(my.ped,all.alleles){
    for(i in 1:nrow(my.ped)){
        my.ind<-unlist(my.ped[i,])
        if(my.ind["mumid"]!=0){
            all.alleles[i,1]<-sample(all.alleles[my.ind["mumid"],],1)
        }
        if(my.ind["dadid"]!=0){
            all.alleles[i,2]<-sample(all.alleles[my.ind["dadid"],],1)
        }
        
    }
    all.alleles
}

sim.alleles.pedigree<-function(my.ped, which.set){
    all.alleles<-matrix(ncol=2,nrow=nrow(my.ped))
    if(length(which.set)==2) all.alleles[which.set,]<-rbind(par.1,par.2)
    if(length(which.set)==1) all.alleles[which.set,]<-rbind(par.1)
    
    sim.alleles<-generate.alleles(my.ped,all.alleles)
    
    missing.genos<-apply(all.alleles,1,function(my.ind){all(is.na(my.ind))})
    all.alleles[missing.genos,]<-sim.alleles[missing.genos,]
    
    
    plot.ped.allele(my.ped,which.allele=1,allele.cols =all.alleles[,1])
    plot.ped.allele(my.ped,which.allele=2,allele.cols=all.alleles[,2])
    plot.ped.allele(my.ped,which.allele=0,allele.cols =rep("black",nrow(all.alleles)))
}

## currently colours thhe alleles of up to 2 inds in pedigree
allele.cols<-brewer.pal(4,"Dark2")
par.1<-allele.cols[1:2]
par.2<-allele.cols[3:4]

############full sib pedigree
sib.ped<-data.frame(id=1:4)    #c("brother","sister","mum","dad"))
sib.ped$dadid<-c(3,3,0,0)    #c("dad","dad",NA,NA)
sib.ped$mumid<- c(4,4,0,0)
sib.ped$ped<-rep(100,4)
sib.ped$sex<-c(1,2,1,2)

############### inbred kid of full sib pedigree
inbred.sib.ped<-rbind(sib.ped,c(5,1,2,100,1))

############### 1/2 sib pedigree

half.sib.ped<-data.frame(id=1:5)  
half.sib.ped$dadid<-c(3,5,0,0,0)   
half.sib.ped$mumid<- c(4,4,0,0,0)
half.sib.ped$ped<-rep(100,5)
half.sib.ped$sex<-c(1,2,1,2,1)


###############inbred kid f 1/2 sib pedigree 
inbred.half.sib.ped<-rbind(half.sib.ped,c(6,1,2,100,1))

############ full cousins

cousin.ped<-data.frame(id=1:8)   
cousin.ped$id<-c(1,2,3,4,5,6,7,8) #c("gp","gm","m1","m2","d1,d2,c1,c2))
cousin.ped$dadid<-c(0,0,2,2,0,0,5,6)    #c("dad","dad",NA,NA)
cousin.ped$mumid<- c(0,0,1,1,0,0,3,4)
cousin.ped$ped<-rep(100,8)
cousin.ped$sex<-c(2,1,2,2,1,1,1,2)

############kid of inbred full cousins
inbred.cousin.ped<-rbind(cousin.ped, c(9,7,8,100,1))
    
run.ped.sim<-function(ped.type, inbreed){
    if(!as.logical(inbreed)){
        if(ped.type=="Full first cousins") sim.alleles.pedigree(cousin.ped,c(1,2))
        if(ped.type=="Full sib") sim.alleles.pedigree(sib.ped,c(3,4))
        if(ped.type=="Half sib")  sim.alleles.pedigree(half.sib.ped,4)
    }
    if(as.logical(inbreed)){
        if(ped.type=="Full first cousins") sim.alleles.pedigree(inbred.cousin.ped,c(1,2))
        if(ped.type=="Full sib") sim.alleles.pedigree(inbred.sib.ped,c(3,4))
        if(ped.type=="Half sib")  sim.alleles.pedigree(inbred.half.sib.ped,4)        
    }

            # if(ped.type=="Half sib")  sim.alleles.pedigree(inbred.sib.ped,c(3,4))
    
}

#user interface
ui <- pageWithSidebar( 
    
    headerPanel = headerPanel("Pedigree Identity by Descent."),
    
    sidebarPanel(
                selectInput("Pedigree_type", "Pedigree_type:",
                        choices = c("Full sib","Half sib","Full first cousins")),
        checkboxInput("Inbreed", "inbred child", FALSE),
        actionButton("goButton", "GO"),
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
        return(c(input$Pedigree_type,input$Inbreed))
    })
    
    
    output$freq <- renderPlot({
        fam_type <- rand()
        run.ped.sim(fam_type[1],fam_type[2])
    })
    

}

# Run the application 
shinyApp(ui = ui, server = server)