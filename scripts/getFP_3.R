# Rscript scripts/timed_getFP.R /data/stampr_read_table-cage1-mA.csv NULL "1" 0.5 /output/timed-test.csv
# to test locally using docker r environment
#
#docker run -it --rm \
#-v /Users/ruprec01/Documents/Faith_lab/Git/bc_seq:/data \
#amplicon-ez /bin/bash
#source activate r_env
#Rscript /data/scripts/getFP_2.R /data/test/total_ST4_readstable.csv NULL "1" 0.03 /data/test/FP_out

##minweight of 0.03 works well
getFP <- function(ReadsTableName, CFUtable, WhereAreReferences, minweight, output_dir) {
  library(EnvStats)
  ReadsTable <- read.csv(ReadsTableName, row.names = 1)
  
  if(!is.null(CFUtable)) {
    CFUtable <- read.csv(CFUtable)
  }
  
  ReferenceVector <- rowSums(cbind(ReadsTable[,WhereAreReferences]))

 
  
  TestNames <- colnames(ReadsTable)
  TableWithoutNoise <- data.frame(row.names = rownames(ReadsTable))
  
  
  
  ResiliencyIndices <- function(samplename, plots = TRUE, ResiliencyLimit = 20000, FractionSD = 20){
    print('analyzing sample')
    print(samplename)

    ###Defining Values###
    invec <- as.numeric(ReferenceVector)
    outvec <- as.numeric(ReadsTable[,which(TestNames == samplename)])
    
    TotalReads <- sum(outvec)
    outvecfreq <- outvec/sum(outvec)
    outvec[which(outvecfreq < 1E-7)] <- 0
    
    if(!is.null(CFUtable)) {
      cfu <<- CFUtable[CFUtable[,1]==samplename,2]
    } else {cfu <- 1E20}
    
    cfu <- as.numeric(cfu)
    if(length(cfu) == 0) {cfu <- 1E20}


    
    ##smoother between 1 and 1.5
    ratio <- sum(outvec > 1) / sum(outvec == 1)
    factor <- max(0, 1.5-ratio)
    removeones <- min(sum(outvec == 1), round(sum(outvec == 1) * factor))

    if ( sum(outvec > 1) >  1.5*sum(outvec == 1) )
    {times <- min(ceiling(cfu), sum(outvec > 0))
    } else  {times <- min(ceiling(cfu), sum(outvec > 0) - removeones)}

    ###Defining greatest frequency change###
    greatestdif <- NULL
    
    #### THIS CREATES ARTIFICIAL VALUES NOT SURE WHY THIS IS WANTED?!
    if(sum(outvec) < 100) {outvec[1:2] <- c(10,100)}
    
    outveccopy <- outvec
    outveccopy[outveccopy < 2] <- 0
    
    outvecsorted <- as.numeric(sort(outveccopy, decreasing = TRUE))
    plotdif <- c(outvecsorted[2:length(outvecsorted)], tail(outvecsorted, 1))
    plotsub <- na.omit(log(outvecsorted) - log(plotdif))
    plotsub <- plotsub[1:length(plotsub)-1]
    
   
    
    ###Second Noise Correction###
    if(sum(outvec == 1) > 1.2* sum(outvec > 1)) {
      outvec <- outvec - 1
      greatestdif <- which(plotsub==max(plotsub))
      outvec[outvec < 0] <- 0
      }

    if(max(plotsub) > 2.302585) {greatestdif <- which(plotsub == max(plotsub))}
    
    ###Adjusting times if there are very few nonzero barcodes###
    if (sum(outvec!=0) < 0.1*sum(invec !=0)) {times <- min(ceiling(cfu), sum(outvec > 1))}
      
    ###Starting Resiliency Plot###
    bind <- cbind(invec, outvec)
    bindsorted <- bind[order(bind[,2]),]
    x <- 0
    outvecsorted <- sort(outvec, decreasing = TRUE)
    
    minusone <- function(dims) {
      input <- bindsorted[,1][1:dims]
      out <- bindsorted[,2][1:dims]
      inputprop <- na.omit(input/sum(input))
      outprop <- na.omit(out / sum(out))
      num <- ( outprop - inputprop ) ^ 2
      den <- inputprop*(1-inputprop)
      sigma <- num / den
      sigma <- sigma[sigma!=Inf]
      F <- mean(na.omit(sigma))
      max( 1 / (F - 1/sum(invec) - 1/sum(outvec)), 1/F)
    }
    
    if(times == 0) {times <- 1}
    timestorep <- rev(seq(from = length(invec) - times + 1, to = length(invec), by = 1))
    
    if(times > ResiliencyLimit) {x <- rep(1, times)} else 
    {x <- sapply(timestorep, minusone)}
    
    
    scanformin <- function(p) {
      start <- p
      findmin <- function() {
        where <- which(x == start)[1]
        #newlocation <- abs(rnorm(1, mean = where, sd = length(na.omit(x))/10))
        newlocation <- abs(rnorm(1, mean = where, sd = sum(outvec>1)/FractionSD))
        if(newlocation > length(na.omit(x))) {newlocation <- where}
        if(newlocation < 1) {newlocation <- where}
        newstart <- x[round(newlocation)]
        if(newstart < start) {start <- newstart}
        start<<-start
      }
      startvector <- replicate(20000, findmin())
      decision <- which(x==start)[1]
    }
    
    q <- seq(1, length(na.omit(x)), length.out = log(length(na.omit(x))))
    p<-x[q]
    
    if(mean(na.omit(x)) == 1) {guesses <- length(na.omit(x))} else
    {guesses <- sapply(p, scanformin)}
    guesses <<- guesses
    
    
    if(length(guesses) == 0) {guesses <- 1}

    ###Defining breaks###
    guessesuniquesorted <- sort(unique(c(guesses)))
    guessesuniquesorted <- c(guessesuniquesorted, length(na.omit(x)))
    
    secondlastguess <- guessesuniquesorted[length(guessesuniquesorted) -1]
    if (rev(sort(outvec))[secondlastguess] - rev(sort(outvec))[max(guessesuniquesorted)] < 4)  {guessesuniquesorted <- guessesuniquesorted[1:length(guessesuniquesorted)-1]}
    
    
     xdif <- c(log(x[1]), log(x))
     xdif <- xdif[1:length(xdif)-1]
     xsub <- (log(x) - xdif)
     greatestdifx <- ifelse(max(xsub) > 2.302585, which(xsub == max(xsub))-1, max(guessesuniquesorted)) 
    
 
    if(is.null(greatestdif)) {greatestdif <- length(na.omit(x))}
    guessesuniquesorted <- sort(unique(c(guessesuniquesorted, greatestdif, greatestdifx)))
    guessesuniquesorted <- guessesuniquesorted[guessesuniquesorted!=0]
    guessesuniquesortedstaggered <- c(-10000000, guessesuniquesorted)
    difference <- c(guessesuniquesorted, 10000000) - guessesuniquesortedstaggered
   
    
     if (min(difference) == 1 & cfu > 2) {
      breaktoremove <- max(which(difference == 1))-1
      guessesuniquesorted <- guessesuniquesorted[-breaktoremove]
    }
    guessesuniquesorted[guessesuniquesorted > max(guessesuniquesorted)-.01*length(na.omit(x))] <- max(guessesuniquesorted)
    guessesuniquesorted <- unique(guessesuniquesorted)
    guessesuniquesorted <<- guessesuniquesorted
    
    ###Creating Indices Table###
    accountsfor <- function(t) {
      outvecsorted <- sort(outvec)
      topnumbers <- tail(outvecsorted, t)
      sum(topnumbers) / sum(outvec)
    }
    
    
    fractionaccounted <- sapply(guessesuniquesorted, accountsfor)
    staggered <- c(0, fractionaccounted)[1:length(fractionaccounted)]
    subtracted <- fractionaccounted - staggered 
    subtracted <- subtracted/sum(subtracted)

    indices <- data.frame(guessesuniquesorted, subtracted)
    colnames(indices) <- c("Number of barcodes", "Fraction of reads")
    
    
    weights <- log(indices[,2])
    values <- (indices[,1])
    weightsforsubtraction <- c(weights[2:length(values)], 0)
    weightsdif <- (weights - weightsforsubtraction)[1:length(weights)-1]
    if(length(weightsdif)==0) {weightsdif <- values}
    noisestart <- indices[,1][which(weightsdif == max(weightsdif))]
    noisestartcopy <- noisestart
    
    #if there are no populations that are below the minimum weight threshold, set noise to be the end of the last detected population
    if (min(indices[,2]) > minweight) {noisestart <- max(indices[,1])}
    #if the sum of the weights of the population after the start of noise is greater than the miniumum weight threshold, set the start to be where the the rest of the reads after are under the minimum weight
    locationofminweightcutoff <- min(which(cumsum(indices[,2])>(1-minweight)))
    if (sum(indices[,2][which(indices[,1] > noisestart)]) > minweight) {noisestart <- indices[,1][locationofminweightcutoff]}
    ##Removing least abundant barcode if it is too nonabundant compared to the second least abundant barcode
    if (noisestart == 2) {
    LeastAbundantBarcode <- sort(outvec, decreasing=TRUE)[noisestart]
    SecondLeastAbundantBarcode <- sort(outvec, decreasing=TRUE)[noisestart-1]
    if(LeastAbundantBarcode/SecondLeastAbundantBarcode < .01) {noisestart <- noisestart - 1}
    }
    
    numberofzeros <- length(invec) - noisestart
    outvecwithoutnoise <- outvec
    outvecwithoutnoise[order(outvecwithoutnoise)][1:numberofzeros] <- 0
    
    ##Removing least abundant barcode if it is too nonabundant compared to the second least abundant barcode
    as.numeric(head(sort(outvecwithoutnoise, decreasing = TRUE)))
    
    FirstResample <- as.numeric(rmultinom(1, sum(outvecwithoutnoise), ReferenceVector/sum(ReferenceVector)))
    
    GetNewBotTable <- function(n) {
      #vec <- as.numeric(rmvhyper(1, FirstResample, n))
      vec <- rmultinom(5, n, FirstResample/sum(FirstResample))
      vec[vec!=0] <- 1
      mean(colSums(vec))
    }
    
    steps1 <- round(seq(from = 1, to = sum(invec != 0), length.out = 100))
    steps2 <- round(seq(from = sum(invec != 0), to = sum(invec != 0)*20, length.out = 100))
    steps <- c(steps1, steps2)
    y <- as.numeric(sapply(steps, GetNewBotTable))

    
    interpol <- approx(x = steps, y = y, n = length(invec))
    xvals <- as.numeric(unlist(interpol[1]))
    yvals <- as.numeric(unlist(interpol[2]))
    
    dfxy2 <- floor(sortedXyData(x = xvals, y = yvals))
    Ns <- NLSstClosestX(dfxy2, noisestart)
    
    inputprop <- na.omit(ReferenceVector/sum(ReferenceVector))
    outprop <- na.omit(outvec / sum(outvec))
    num <- ( outprop - inputprop ) ^ 2
    den <- inputprop*(1-inputprop)
    sigma <- num / den
    sigma <- sigma[sigma!=Inf]
    F <- mean(na.omit(sigma))
    Nb <- 1 / (F - 1/sum(invec) - 1/sum(outvec)) 
    AverageFreq <- 1/geoMean(outvecwithoutnoise[outvecwithoutnoise>0]/sum(outvecwithoutnoise[outvecwithoutnoise>0]))
    
    
    ReadsAtWeightCutoff <- (1-minweight)*sum(outvec)
    MinCutoff <- outvecsorted[min(which(cumsum(outvecsorted[outvecsorted > 0]) > ReadsAtWeightCutoff))]
    NBarcodesAtMinweight <- sum(outvec > MinCutoff) + 1  
    Ns_MinCutoff <- NLSstClosestX(dfxy2, NBarcodesAtMinweight)
    
    if(plots){
  # Format the statistics for inclusion in the plot
  stats_text <- paste0(
    "Nb = ", round(Nb, 2), "\n",
    "Ns = ", round(Ns, 2), "\n", 
    "AverageFreq = ", round(AverageFreq, 2), "\n",
    "Ns at minweight = ", round(Ns_MinCutoff, 2)
  )
  
  # Create SVG file
  png(file.path(output_dir, paste0(samplename, "_plot.png")), width = 800, height = 800, res = 300) # Increased height for annotations
  
  # Setup layout with space for annotations
  layout(matrix(c(1, 2, 3), nrow = 3, ncol = 1, byrow = TRUE), heights = c(5, 5, 2))
  
  # Define the cutoff positions function
  converttocutoffs <- function(p) {
    outvecsorted <- sort(outvec, decreasing = TRUE)
    cutoffvalue <- outvecsorted[p]
    if(cutoffvalue == 0) {cutoffvalue <- 1}    
    cutoffvalue/sum(outvec) 
  }
  cutoffpositions <- unlist(lapply(guessesuniquesorted, converttocutoffs))
  cutoffpositions <<- cutoffpositions
  
  # First plot - Resiliency
  par(mar = c(4, 4, 3, 2)) # Adjust margins
  plot(1:times, x, log = "y", ylim = c(.01, 2E6), xlim = c(1, length(as.numeric(na.omit(x)))), 
       main = "Resiliency", ylab = "Nb", xlab = "Iteration")
  
  # Add vertical lines with labels
  for (i in seq_along(guessesuniquesorted)) {
    abline(v = guessesuniquesorted[i], col = rainbow(length(guessesuniquesorted))[i], lwd = 1.5)
    # Add label for each line
    text(guessesuniquesorted[i], max(x, na.rm=TRUE) * 0.8, 
         labels = paste0("G", i, ": ", guessesuniquesorted[i]), 
         srt = 90, adj = c(1.1, 0.5), cex = 0.7, 
         col = rainbow(length(guessesuniquesorted))[i])
  }
  
  # Second plot - Frequency
  par(mar = c(4, 4, 3, 2)) # Adjust margins
  plot(outvec/sum(outvec), log = "y", ylim = c(1E-8, 1), 
       main = samplename, ylab = "Frequency", xlab = "Barcode")
  
  # Add horizontal lines with labels
  for (i in seq_along(cutoffpositions)) {
    abline(h = cutoffpositions[i], col = rainbow(length(cutoffpositions))[i], lwd = 1.5)
    # Add label for each line
    text(length(outvec) * 0.9, cutoffpositions[i] * 1.2, 
         labels = paste0("C", i, ": ", round(cutoffpositions[i], 7)), 
         adj = c(1, 0.5), cex = 0.7, 
         col = rainbow(length(cutoffpositions))[i])
  }
  
  # Third panel - Statistics text
  par(mar = c(0, 1, 0, 1)) # Minimal margins for text panel
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(0, 0, stats_text, cex = 0.9)
  
  # Print indices in a separate table
  textplot_indices <- function(df) {
    # Create a simple text representation of the dataframe
    colnames(df) <- c("Num Barcodes", "Frac Reads")
    
    # Format the numbers to be more compact
    df_formatted <- data.frame(
      "Num Barcodes" = format(df[,1], scientific = FALSE),
      "Frac Reads" = format(df[,2], digits = 3, scientific = FALSE)
    )
    
    text_table <- capture.output(print(df_formatted, row.names = FALSE))
    text_table <- paste(text_table, collapse = "\n")
    
    # Print the table to the current device
    mtext(text_table, side = 1, line = -7, at = 0, adj = 0, cex = 0.8)
  }
  
  # Add indices table below stats
  if(nrow(indices) <= 10) { # Only show if not too large
    textplot_indices(indices)
  }
  
  # Close device
  dev.off()
  
  # Still print to console
  print(indices)
  print(paste("Nb=", Nb, sep = ""))
  print(paste("Ns=", Ns, sep = ""))
  print(paste("AverageFreq=", AverageFreq, sep = ""))
  print(paste("Ns at minweight=", Ns_MinCutoff, sep= ""))
  
  # Assign global variables as before
  Nb <<- Nb
  Ns <<- Ns
  AverageFreq <<- AverageFreq
  indices <<- indices
  outvec <<- outvec
  invec <<- invec
  dfxy2 <<- dfxy2
  FirstResample <<- FirstResample
  noisestart <<- noisestart
  x <<- x
  greatestdif <<- greatestdif
  cfu <<- cfu
  times <<- times
  originaloutvec <<- ReadsTable[,which(TestNames == samplename)]
  MinCutoff <<- MinCutoff
}
    
    if(plots == TRUE) {
      print(samplename)
      print(indices)
      print(Ns)
      
      numberofzeros <- length(invec) - noisestart
      outvec[order(outvec)][1:numberofzeros] <- 0
      
      TableWithoutNoise <- data.frame(TableWithoutNoise, outvec)
      colnames(TableWithoutNoise)[length(colnames(TableWithoutNoise))] <- samplename
      TableWithoutNoise <<- TableWithoutNoise
      outvec <<- outvec
      invec <<- ReferenceVector
      if(cfu == 1E20) {cfu <- 0}
      
      den <- density(dfxy2$y)
      denapprox <- approx(den, n=length(invec))
      deny <- denapprox$y
      scalefac <- 1/max(deny)
      uncertainty <- (deny*scalefac)[noisestart]
      
      c(TotalReads, noisestart, Ns_MinCutoff, Nb, Ns, AverageFreq, cfu, log10(Ns), log10(cfu), cfu/Ns)
    }
  }
  

  SampleNames <- colnames(ReadsTable)[-WhereAreReferences]
  TableOfEstimates <- t(sapply(SampleNames, ResiliencyIndices))
  print(TableOfEstimates)
  colnames(TableOfEstimates) <- c("TotalReads", "Number of barcodes", "Ns_MinCutoff", "Nb", "Ns","AverageFrequency", "CFU", "Log10Ns", "Log10CFU", "CFU/Ns")
  print(TableOfEstimates)
  print('hello')
  write.csv(TableOfEstimates, paste0(output_dir, '/TableOfEstimates.csv'))
  write.csv(TableWithoutNoise, paste0(output_dir, "/FrequenciesWithoutNoise.csv"))
  TableOfEstimates <<- TableOfEstimates
  FrequenciesWithoutNoise <<- TableWithoutNoise
  ReadsTable <<- ReadsTable
  CFUtable <<- CFUtable
  ReferenceVector <<- ReferenceVector
  TestNames <<- TestNames
  minweight <<- minweight
  ReadsTableName <<- ReadsTableName
  WhereAreReferences <<- WhereAreReferences
}
cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) != 5) {
  stop("Usage: Rscript getFP_12112023.R <ReadsTable> <CFUtable> <WhereAreReferences> <minweight> <output_dir>")
}

ReadsTable <- cmd_args[1]
CFUtable <- cmd_args[2]
if (CFUtable == "NULL") {
  CFUtable <- NULL
}
WhereAreReferences <- as.integer(unlist(strsplit(cmd_args[3], ",")))
minweight <- as.numeric(cmd_args[4])
output_dir <- cmd_args[5]

# Print parsed arguments for debugging
print("Parsed arguments:")
print(paste("ReadsTable:", ReadsTable))
print(paste("CFUtable:", CFUtable))
print(paste("WhereAreReferences:", WhereAreReferences))
print(paste("minweight:", minweight))
print(paste("output_dir:", output_dir))

# Call the getFP function with the provided arguments
getFP(ReadsTable, CFUtable, WhereAreReferences, minweight, output_dir)