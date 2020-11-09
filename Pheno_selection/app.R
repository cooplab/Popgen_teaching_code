#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

##Paste this in first
par.off.corr<-function(L=20, environ.var,Num_inds=1000,print.slope=FALSE,selection=FALSE,sel.cutoff){
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
    if(selection){
        sel.cutoff<-quantile(c(dad.pheno,mum.pheno),p=1-sel.cutoff)
        plot(parental.midpoint,child.pheno,xlab="Parental midpoint",ylab="Child's phenotype",cex=1.5,cex.axis=1.5,cex.main=1.5,cex.lab=1.5,
             main=paste("Mid-par offspring regression slope ",format( my.slope,digit=2)),col=ifelse(parental.midpoint>sel.cutoff,"red","grey"),asp = 1)
    }else{
        plot(parental.midpoint,child.pheno,xlab="Parental midpoint",ylab="Child's phenotype",cex=1.5,cex.axis=1.5,cex.main=1.5,cex.lab=1.5,
             main=paste("Mid-par offspring regression slope ",format( my.slope,digit=2)),asp = 1)
    }
    ## plot the regression in red
    abline(h=0,col="grey",lwd=2)
    abline(v=0,col="grey",lwd=2)
    abline(lm(child.pheno~parental.midpoint),col="blue",lwd=2)
    
    
    
    
    if(selection){
        sel.child.mean<-mean(child.pheno[parental.midpoint>sel.cutoff])
        sel.par.mean<-mean(parental.midpoint[parental.midpoint>sel.cutoff]);
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


one.gen.sel<-function(L=1000,environ.var,sel,plot.geno=FALSE,add.arrows=FALSE){
    ##Quantitative genetics sims
    allele.freq<-0.5   ###each locus is assumed to have the same allele frequencies. This is just to simplify the coding, in reality these results work when each locus has its own frequency (and the coding wouldn't be too much harder). 
    
    
    Num_inds=10000
    
    ##MAKE A MUM
    ## For each mother, at each locus we draw an allele (either 0 or 1) from the population allele frequency. 
    ##We do this twice for each mother two represent the two haplotypes in the mother 
    mum.hap.1<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    mum.hap.2<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    ##type mum.hap.1[,1] to see the 1st mothers 1st haplotype
    
    ##Each mothers genotype at each locus is either 0,1,2
    mum.geno<-mum.hap.1+mum.hap.2
    
    additive.genetic<-colSums(mum.geno)
    mean.genetic<-mean(additive.genetic)
    genetic.var<-sd(additive.genetic)
    
    additive.genetic<-additive.genetic / sd(additive.genetic)
    mum.pheno<- additive.genetic + rnorm(Num_inds,sd=sqrt(environ.var))
    mum.pheno<-mum.pheno-mean(mum.pheno)
    
    
    
    ###FAMILIES
    
    
    ##MAKE A DAD (same code as make a mum
    dad.hap.1<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    dad.hap.2<-replicate(Num_inds, rbinom(L,1,allele.freq) )
    dad.geno<-dad.hap.1+dad.hap.2
    
    
    additive.genetic<-colSums(dad.geno)
    additive.genetic<-additive.genetic / sd(additive.genetic)
    dad.pheno<- additive.genetic + rnorm(Num_inds,sd=sqrt(environ.var))
    dad.pheno<-dad.pheno-mean(dad.pheno)
    
    ### Make a child
    child.geno<-dad.hap.1+mum.hap.1 ##1/2 from mum 1/2 from dad
    
    additive.genetic<-colSums(child.geno)
    additive.genetic<-additive.genetic / sd(additive.genetic)
    child.pheno<- additive.genetic + rnorm(Num_inds,sd=sqrt(environ.var))
    child.pheno<-child.pheno-mean(child.pheno)
    
    
    
    ##Selection of top sel% of individuals
    
    top.sel.per.mums<- mum.pheno>quantile(mum.pheno,p=1-sel) 
    top.sel.per.dads<- dad.pheno>quantile(dad.pheno,p=1-sel)
    
    
    child.geno<-dad.hap.1[,top.sel.per.dads]+mum.hap.1[,top.sel.per.mums] ##1/2 from mum 1/2 from dad
    
    additive.genetic<-(colSums(child.geno)-mean.genetic)
    additive.genetic<-additive.genetic/genetic.var
    child.pheno<- additive.genetic + rnorm(length(child.geno),sd=sqrt(environ.var))
    
    layout(1:3)
    my.lim<-quantile(c(mum.pheno,dad.pheno),p=c(0.01,0.99))
    my.lim[2]<-quantile(child.pheno,p=c(0.99))
    
    hist(c(mum.pheno,dad.pheno),breaks=100,xlim=my.lim,xlab="Phenotype",main=paste("Phenotype distribution before selection"),cex.axis=1.5,cex.lab=1.5,cex.main=1.5); #, Mean=0, VA=1, VE=",environ.var,", Taking top ",round(100*sel),"%",sep="")
    abline(v=0,col="blue",lwd=3)
    
    par.mean<-mean(c(mum.pheno[top.sel.per.mums],dad.pheno[top.sel.per.dads]))
    hist(c(mum.pheno[top.sel.per.mums],dad.pheno[top.sel.per.dads]),breaks=100,xlim=my.lim,xlab="Phenotype",main=paste("Phenotype distribution after selection, parental mean=",format(par.mean,dig=3)),cex.axis=1.5,cex.lab=1.5,cex.main=1.5); 
    abline(v= par.mean,col="red",lwd=3)
    abline(v=0,col="blue",lwd=3)
    
    
    if(add.arrows){
        arrows(x0=0,x1=par.mean,y0=50,y1=50,col="blue",lwd=2,length=0.1,code=3)
        text(x=par.mean/2, y=70,"S",col="blue",cex=1.5)
    }
    
    hist(child.pheno,xlim=my.lim,breaks=100,xlab="Phenotype",main=paste("Phenotype distribution in the children Mean in children = ",format(mean(child.pheno),dig=3)),cex.axis=1.5,cex.lab=1.5,cex.main=1.5); 
    abline(v=0,col="blue",lwd=3)
    abline(v= mean(child.pheno),col="red",lwd=3)
    
    if(add.arrows){
        arrows(x0=0,x1=mean(child.pheno),y0=500,y1=500,col="blue",lwd=2,length=0.1,code=3)
        text(x=mean(child.pheno)/2, y=800,"R",col="blue",cex=1.5)
    }
    
    ##Mean phenotype after selection
    cat("Selected parental mean",par.mean,"\n")
    ##Mean child phenotype
    cat("Mean in children = ",mean(child.pheno),"\n")
    
    if(plot.geno){
        #	quartz()
        layout(1:2)
        par(mar=c(4,4,1.5,1))
        sel.dad.genosum<-colSums(dad.geno[,top.sel.per.dads])
        rand.dad.genosum<-colSums(dad.geno[,sample(top.sel.per.dads)])
        
        a<-hist(rand.dad.genosum,breaks=20,plot=FALSE)
        b<-hist(sel.dad.genosum,breaks=20,plot=FALSE)
        #recover()
        y.max<-max(c(a$counts,b$counts))
        hist(rand.dad.genosum,breaks=20,  col = rgb ( 1 , 0 , 0 , 0.4 ),xlim=c(min(rand.dad.genosum)-5,max(sel.dad.genosum)+5),,ylim=c(0, y.max),xlab="Num. up alleles", main="Parental generation",cex.lab=1.5)
        hist(sel.dad.genosum,breaks=20,  col =rgb ( 0 , 0 , 1 , 0.4 ),add=TRUE) #
        cat("Number up alleles in selected pars ",mean(colSums(dad.geno[,top.sel.per.dads])))
        legend ( "topleft" , legend = c ( "All individuals" , "Selected parents" ) , pch = 15 , col = c ( rgb ( 1 , 0 , 0 , 0.4 )  , rgb ( 0 , 0 , 1 , 0.4 ) ) , bty = "n" , cex = 1.5 )
        ###hist of kids
        hist(colSums(child.geno),breaks=20,  col =rgb ( 0 , 1 , 0 , 0.4 ),xlim=c(min(rand.dad.genosum)-5,max(sel.dad.genosum)+5),xlab="Num. up alleles", main="Next generation",cex.lab=1.5)
        cat("Number up alleles in kids ",mean(colSums(child.geno)))
        
        legend ( "topleft" , legend = c ( "Children" ) , pch = 15 , col = rgb ( 0 , 1 , 0 , 0.4 ), bty = "n" , cex = 1.5 )
    }
}



#user interface
ui <- pageWithSidebar( 
    
    headerPanel = headerPanel("Phenotypic selection"),
    
    sidebarPanel(
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
         sliderInput("selcutoff","Truncation (% of pop. selected): ",
                     min=5,
                     max=50,
                     value=50),
    
        actionButton("goButton", "GO"),
    ), 
    #mainPanel =  mainPanel(

    mainPanel(type = "tabs",
        tabsetPanel(
            tabPanel("Pheno. distribution", 
                     helpText("
                     The three tabs above show three different views of selection on the same trait.
                     The graph below shows the distribution of phenotypes before and after phenotypic selection. 
                              You get to choose what percentage of the population is selected to form the next generation. 
                              Our VA is set to 1, and you get to vary VE and the number of loci underlying the trait."),
                     plotOutput(outputId = 'pheno')
            ),
            tabPanel("Geno. distribution", 
                     helpText("The distribution of genotypes before and after phenotypic selection in the population. 
                     The truncation selection on the phenotype results in a shift in the genotypes of the population. 
                              The histograms below the number of trait increasing alleles carried by individuals.
                              Note how the children of the selected individual are faithful representation of their parents." ),
                     plotOutput(outputId = 'geno')
            ),
            tabPanel("Parent-Offspring plot",
                     helpText("The covariance of parental mid-point and offspring  phenotype. 
                              The complete set of dots (red and grey) shows the relationship in the absence of selection. 
                              The red dots show the relationship for those parental individuals selected to reduce under truncation selection"),
                     
                     plotOutput(outputId = 'par_off')  
            )
        )
    )
)
# Define server logic 
#back end code and response to user input
server <- function(input, output){
    
    #each time user hits "go"
    rand <- eventReactive(input$goButton, {
        #parameters
        
        #package data for plotting
        return(input)
    })
    
    output$par_off <- renderPlot({
        my.params<- rand()
        par.off.corr(L=my.params$L,environ.var=my.params$VE,Num_inds=1000,print.slope=FALSE,
                               selection=TRUE,sel.cutoff=my.params$selcutoff/100)

    })    
        
    output$pheno <- renderPlot({
        my.params<- rand()
        one.gen.sel(L=my.params$L,environ.var=my.params$VE,
                              sel=my.params$selcutoff/100,plot.geno=FALSE,add.arrows=TRUE)
    })
    
    output$geno <- renderPlot({
        my.params<- rand()
        one.gen.sel(L=my.params$L,environ.var=my.params$VE,
                    sel=my.params$selcutoff/100,plot.geno=TRUE,add.arrows=FALSE)
    })
}



# Run the application 
shinyApp(ui = ui, server = server)
