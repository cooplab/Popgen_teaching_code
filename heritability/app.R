#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)


IBD.coeffs<-data.frame(IBD0=0,IBD=0,IBD2=0)
IBD.coeffs["Identical Twins",]<-c(0.,0.,1)
IBD.coeffs["Parent & Child",]<-c(0.0,1.0,0.0)
IBD.coeffs["Full Sibs",]<-c(0.25,0.5,0.25)
IBD.coeffs["1/2 Sibs",] <-c(0.5,0.5,0)
IBD.coeffs["1st Cousins",] <-c(0.75,0.25,0.0)

##Paste this in first
par.off.corr<-function(L=20, environ.var,Num_inds=1000,print.slope=FALSE,sel.cutoff=FALSE){
    ##Quantitative genetics sims
    allele.freq<-0.5   ###each locus is assumed to have the same allele frequencies. This is just to simplify the coding, in reality these results work when each locus has its own frequency (and the coding wouldn't be too much harder). 
    
    
    ##MAKE A MUM
    ## For each mother, at each locus we draw an allele (either 0 or 1) from the population allele frequency. 
    ##We do this twice for each mother two represent the two haplotypes in the mother 
    mum.hap.1<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    mum.hap.2<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    ##type mum.hap.1[,1] to see the 1st mothers 1st haplotype
    
    ##Each mothers genotype at each locus is either 0,1,2
    mum.geno<-mum.hap.1+mum.hap.2
    
    additive.genetic<-colSums(mum.geno)
    genetic.sd<-sd(additive.genetic)
    mean.genetic<-mean(additive.genetic)
    
    additive.genetic<-additive.genetic / sd(additive.genetic)
    mum.pheno<- additive.genetic + rnorm(Num_inds,sd=sqrt(environ.var))
    mum.pheno<-mum.pheno-mean(mum.pheno)
    
    
    
    ##MAKE A DAD (same code as make a mum
    dad.hap.1<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    dad.hap.2<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    dad.geno<-dad.hap.1+dad.hap.2
    
    
    additive.genetic<-colSums(dad.geno)
    additive.genetic<-additive.genetic / sd(additive.genetic)
    dad.pheno<- additive.genetic + rnorm(Num_inds,sd=sqrt(environ.var))
    dad.pheno<-dad.pheno-mean(dad.pheno)
    
    ### Make a child
    child.geno<-dad.hap.1+mum.hap.1 ##1 haplotype from mum 1 haplotype from dad
    
    additive.genetic<-colSums(child.geno)
    additive.genetic<-additive.genetic / sd(additive.genetic)
    child.pheno<- additive.genetic + rnorm(Num_inds,sd=sqrt(environ.var))
    child.pheno<-child.pheno-mean(child.pheno)
    
    
    
    ##Calculate midpoints, linear model and plots
    
    parental.midpoint<-(mum.pheno+dad.pheno)/2 ##avg. parents
    
    lm(child.pheno~parental.midpoint) ##linear model between child and mid point
    ##the slope of this line is the narrow sense heritability.
    
    my.slope<-lm(child.pheno~parental.midpoint)$coeff[2]
    # plot parental midpoint against offsprings phenotype.
    #layout(1) ###done in case this is run after the code with 3 plots
    if(sel.cutoff){
        plot(parental.midpoint,child.pheno,xlab="Parental midpoint",ylab="Child's phenotype",cex=1.5,cex.axis=1.5,cex.main=1.5,cex.lab=1.5,
             main=paste("Mid-par offspring regression slope ",format( my.slope,digit=2)),col=ifelse(parental.midpoint>1,"red","grey"),asp = 1)
    }else{
        plot(parental.midpoint,child.pheno,xlab="Parental midpoint",ylab="Child's phenotype",cex=1.5,cex.axis=1.5,cex.main=1.5,cex.lab=1.5,
             main=paste("Mid-par offspring regression slope ",format( my.slope,digit=2)),asp = 1)
    }
    ## plot the regression in red
    abline(h=0,col="grey",lwd=2)
    abline(v=0,col="grey",lwd=2)
    abline(lm(child.pheno~parental.midpoint),col="blue",lwd=2)
    
    

    
    if(sel.cutoff){
        sel.child.mean<-mean(child.pheno[parental.midpoint>1])
        sel.par.mean<-mean(parental.midpoint[parental.midpoint>1]);
        points(sel.par.mean,0,col="blue",pch=19,cex=1.2)
        pred.sel.child<-my.slope*sel.par.mean
        points(0,pred.sel.child,col="blue",pch=19,cex=1.2)
        lines(c(0,sel.par.mean), rep(pred.sel.child,2),col="blue",lwd=2)
        lines(rep(sel.par.mean,2), c(0,pred.sel.child),col="blue",lwd=2)
        arrows(x0=0,x1=0,y0=0,y1=pred.sel.child,col="blue",lwd=2,length=0.1,code=3)
        arrows(x0=0,x1=sel.par.mean,y0=0,y1=0,col="blue",lwd=2,length=0.1,code=3)
        text(x=sel.par.mean/2,y=min(child.pheno)*0.06,"S",cex=1.5,col="blue")
        text(x=min(parental.midpoint)*0.06,y=pred.sel.child/2,"R",cex=1.5,col="blue")
    }
    
    if(print.slope) text(x=min(parental.midpoint)*.8,y=max(child.pheno)*.9,label=paste("slope= ",format(my.slope,digit=3)),col="red",lwd=4,cex=1.5)
    
    abline(0,1,col="red",lwd=3,lty=2)
    #	recover()
}

##Paste this in first
genetic.covar<-function(L=20, environ.var,Num_inds=5000,print.slope=FALSE,sel.cutoff=FALSE,ibd.prob,relly.type=""){
    ##Quantitative genetics sims
    allele.freq<-0.5   ###each locus is assumed to have the same allele frequencies. This is just to simplify the coding, in reality these results work when each locus has its own frequency (and the coding wouldn't be too much harder). 
    stopifnot(sum(ibd.prob)==1)
    
    
    ##MAKE A IND 1ss
    ## For each ind, at each locus we draw an allele (either 0 or 1) from the population allele frequency. 
    ##We do this twice for each ind two represent the two haplotypes in the mother 
    ind.hap.1<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    ind.hap.2<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    ##type mum.hap.1[,1] to see the 1st mothers 1st haplotype
    
    ##Each mothers genotype at each locus is either 0,1,2
    ind.geno<-ind.hap.1+ind.hap.2
    
    additive.genetic<-colSums(ind.geno)
    genetic.sd<-sd(additive.genetic)
    mean.genetic<-mean(additive.genetic)
    
    additive.genetic<-additive.genetic / sd(additive.genetic)
    ind.pheno<- additive.genetic + rnorm(Num_inds,sd=sqrt(environ.var))
    ind.pheno<-ind.pheno-mean(ind.pheno)
    
    ##MAKE A IND 2's
    ###routine to generate genotypes for set of 2nd individual based on our genotypes for first
    other.ind.geno<-sapply(1:Num_inds,function(ind){
        sapply(1:L,function(snp){
            num.ibd<-sample(0:2,1,prob=ibd.prob)
            if(num.ibd==0){my.geno<-sum(rbinom(2,1,allele.freq))} 
            if(num.ibd==1){my.geno<-ind.hap.1[snp,ind]+rbinom(1,1,allele.freq)} 
            if(num.ibd==2){my.geno<-ind.geno[snp,ind]}
            return(my.geno)
        })
    })
    
    other.ind.additive.genetic<-colSums(other.ind.geno)
    genetic.sd<-sd(other.ind.additive.genetic)
    mean.genetic<-mean(other.ind.additive.genetic)	
    other.ind.additive.genetic<-other.ind.additive.genetic / sd(other.ind.additive.genetic)
    other.pheno<- other.ind.additive.genetic + rnorm(Num_inds,sd=sqrt(environ.var))
    other.pheno<-other.pheno-mean(other.pheno)
    my.cov<-cov(ind.pheno,other.pheno);
    my.var<-var(c(ind.pheno,other.pheno))
    plot(ind.pheno,other.pheno,xlab="Ind 1's phenotype",ylab="Ind 2's phenotype",
         cex=1.5,cex.axis=1.5,cex.main=1.5,cex.lab=1.5,
         main=paste(relly.type,"Cov=",format(my.cov,digit=3)," Var=",format(my.var,digit=3),sep=" "),col=adjustcolor("black",0.4),asp = 1) #,"L =",L,"VE=",environ.var,"VA=1",sep=", "))
    abline(h=0,col="grey",lwd=2)
    abline(v=0,col="grey",lwd=2)
    abline(lm(other.pheno~ind.pheno),col="blue",lwd=2)
    abline(0,1,col="red",lwd=3,lty=2)
    #	textbox(x=min(ind.pheno)*.5,y=max(other.pheno)*.9,textlist=c("Cov= ",format(my.cov,digit=3)),col="red"))
    #legend(x="topleft",legend=paste("Cov= ",format(my.cov,digit=2)),col="red",lwd=4,cex=1.5,bg="white",pch=NA,lty=NA)
    #text(x=min(ind.pheno)*.5,y=max(other.pheno)*.9,label=paste("Cov= ",format(my.cov,digit=3)),col="red",lwd=4,cex=1.5)
    cat("pheno. covariance=",my.cov,"\n")
    cat("Expected covar=",sum(ibd.prob*c(0,0.5,1)),"\n")
    VA<-my.cov/sum(ibd.prob*c(0,0.5,1))
    cat("VA= ", VA,"h2= ",VA/var(c(ind.pheno,other.pheno)),"\n")
}


#user interface
ui <- pageWithSidebar( 
    
    headerPanel = headerPanel("Heritability and selection simulations"),
    
    sidebarPanel(
        # Sidebar with a slider input for number of bins 
        selectInput("Relationship", "Relationship:",
                    choices = c("Mid-parent Offspring regression",
                                "Parent & Child",
                                "Identical Twins",
                                "Full Sibs",
                                "1/2 Sibs",
                                "1st Cousins")),
        sliderInput("L",
                    "Num. loci:",
                    min = 2,
                    max = 100,
                    value = 30),
        sliderInput("VE",
                    "Environmental variance",
                    min = 0.01,
                    max = 5,
                    value = 1),
        # conditionalPanel(condition = "input.Relationship == \"Mid-parent Offspring regression\"",
        #         checkboxInput("Selection", "Selection", FALSE),
        # ),
        # conditionalPanel(condition = "input.Selection == true",
        # sliderInput("selcutoff","Truncation: ",
        #             min=-2.1,
        #             max=3.1,
        #             value=1)
        # ),
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
        if(my.params$Relationship=="Mid-parent Offspring regression"){ 
  #          if(!my.params$Selection){ 
                par.off.corr(L=my.params$L, environ.var=my.params$VE,Num_inds=500)  #,print.slope=TRUE)
  #          }else{
  #              par.off.corr(L=my.params$L, environ.var=my.params$VE,Num_inds=500,sel.cutoff=my.params$selcutoff)  #,print.slope=TRUE)
  #          }
            }else{
            genetic.covar(L=my.params$L, environ.var=my.params$VE,Num_inds=500,ibd.prob=IBD.coeffs[my.params$Relationship,],relly.type=my.params$Relationship)
        }
    
        
    })
    
    
}



# Run the application 
shinyApp(ui = ui, server = server)

