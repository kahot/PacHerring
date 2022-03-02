
plot_StackedRuns_ph <- function(runs=runs, savePlots=FALSE, separatePlots=FALSE, outputName=NULL) {
    
    # avoid notes
    chrom <- NULL
    from <- NULL
    to <- NULL
    
    runs<-runs
    plot_list <- list()
    #select a POPULATION
    for (g in unique(runs$id)){
        print(paste('Current population: ',g))
        teilsatz <- subset(runs,runs$id==g)
        
        chr_order <- c(0:26)
        list_chr=unique(teilsatz$chrom)
        new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))
        
        #select a chromosome
        for (chromosome in new_list_chr){
            
            print(paste('CHR: ',chromosome))
            krom <- subset(teilsatz,chrom==chromosome)
            krom <- krom[order(krom$from),]
            
            #start the order
            yread <- c(); #keeps track of the x space that is used up by segments
            
            # get x axis limits
            minstart <- min(krom$from);
            maxend <- max(krom$to);
            
            # initialise yread
            yread[1] <- minstart - 1;
            ypos <- c(); #holds the y pos of the ith segment
            
            for (r in 1:nrow(krom)){
                read <- krom[r,];
                start <- read$from;
                placed <- FALSE;
                
                # iterate through yread to find the next availible
                # y pos at this x pos (start)
                y <- 1;
                while(!placed){
                    
                    if(yread[y] < start){
                        ypos[r] <- y;
                        yread[y] <- read$to;
                        placed <- TRUE;
                    }
                    
                    # current y pos is used by another segment, increment
                    y <- y + 1;
                    # initialize another y pos if we're at the end of the list
                    if(y > length(yread)){
                        yread[y] <- minstart-1;
                    }
                }
            }
            
            maxy <- length(yread);
            krom$ypos <- ypos;
            utils::head(krom)
            
            #PLOT STACKED RUNS
            p <- ggplot2::ggplot()
            p <- p + ggplot2::geom_segment(data=krom, aes(x = from/(10^6), y = ypos, xend = to/(10^6), yend = ypos),
                                           colour="lightcoral", alpha=1, size=0.75)
            p <- p + xlim(0, max(krom$to/(10^6))+10) + ylim(0,length(yread)+1)
            p <- p + ylab('n Runs') + xlab('Chromosome position (Mbps)')
            p <- p + ggplot2::ggtitle(paste("Group: ",g,'\nChromosome:',chromosome))
            p <- p + theme(plot.title = element_text(hjust = 0.5))
            
            # Save plots by Chromosome
            if(savePlots & separatePlots) {
                if (! is.null(outputName)) { fileNameOutput <- paste(outputName, "Chr", chromosome, r, "Stacked", sep="_")
                } else { fileNameOutput <- paste("Runs_StackedChr", chromosome, r,sep="_") }
                ggsave(filename = paste(fileNameOutput,'.pdf',sep='') , plot = p, device = "pdf")
            } else if (savePlots) { plot_list[[chromosome]] <- p
            } else { print(p) }
        }
        
        # Save plot all Chromosome
        if(savePlots & !separatePlots) {
            if (! is.null(outputName)) { fileNameOutput <- paste(outputName, g ,"StackedAllChr.pdf", sep="_")
            } else { fileNameOutput <- paste("Runs_StackedAllChr",g,".pdf",sep='_') }
            plot_list_final <- gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1)
            ggsave(filename = fileNameOutput , plot = plot_list_final, device = "pdf")
        }
    }
}


