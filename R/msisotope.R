isotope.dist <- function(elements.count) {
  elements <- c( "C", "H", "N", "O", "S" )
  heavy <- c(1.10, 0.015, 0.37, 0.20, 4.21)/100
  names(heavy) <- elements
  light <- 1.00 - heavy
  names(elements.count) <- elements
  single.prob <- as.list( elements )
  names(single.prob) <- elements
  all.prob <- numeric(0)
  for ( e in elements ) {
    count <- elements.count[e]
    v <- seq(0,count)
    l <- light[e]
    h <- heavy[e]
    new.prob <- single.prob[[e]] <- round(choose(count,v)*(l^(count-v))*(h^v),4)
    if (e =="O" | e =="S") { # O and S isotopes are 2 Da more
      new.prob <- rep(0,2*count+1)
      for( i in 1:(count+1)) {
        new.prob[(i-1)*2+1] <- single.prob[[e]][i]
      }
      single.prob[[e]] <- new.prob
    }
    #print(single.prob[[e]])
    if ( length(all.prob) > 0 ) {
      d <- length(all.prob)-length(new.prob)
      if ( d > 0) {
        new.prob <- c(new.prob,rep(0,d))
      } else if ( d< 0 ) {
        all.prob <- c( all.prob,rep(0,-d) )
      }
      all.prob <- round(convolve(all.prob, rev(new.prob),type="o"),4)
    } else {
      all.prob <- single.prob[[e]]
    }
  }
  return(all.prob)
##  return(single.prob)
}

## isotope.dist(c(60,13,13,86,2))
averagine.count <- function(input.mass) {
  averagine.mass <- 111.1254
  elements <- c( "C", "H", "N", "O", "S" )
  averagine.comp <- c( 4.9348, 7.7583, 1.3577, 1.4773, 0.0417 )
  names(averagine.comp) <- elements
  return( round(averagine.comp*(input.mass/averagine.mass)) )
}

## ms2 triggered?
find.ms2.triggered <- function(xfile, yfile, predicted.mz, rt.range) {
  ms1.scanNums <- which( xfile@scantime>=rt.range[1]&xfile@scantime<=rt.range[2] )
  ms2.matrix <- matrix(0, nrow=length(ms1.scanNums), ncol=length(predicted.mz) )
  dimnames(ms2.matrix)[[1]] <- as.character(ms1.scanNums)
  dimnames(ms2.matrix)[[2]] <- as.character(predicted.mz)
  ms2.acq.num.max <- max(xfile@msnAcquisitionNum)
  for (i in 1:dim(ms2.matrix)[1]) {
    ms1.scanNum <- ms1.scanNums[i]
    ms2.acq.num.begin <- xfile@acquisitionNum[ms1.scanNum]+1
    if ( ms1.scanNum < length(xfile@scantime) ) {
      ms2.acq.num.end <- xfile@acquisitionNum[ms1.scanNum+1]-1
    } else {
      ms2.acq.num.end <- ms2.acq.num.max
    }
    ms2.mz <- yfile[yfile[,"num"]>=ms2.acq.num.begin&yfile[,"num"]<=ms2.acq.num.end, "pmz"]
    for (j in 1:dim(ms2.matrix)[2]) {
      ms2.matrix[i,j] <- sum( abs( ms2.mz-predicted.mz[j] )<0.1)
    }
  }
  return(ms2.matrix)
}
