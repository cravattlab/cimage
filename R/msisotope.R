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
##
estimateChromPeakNoise <- function(intensity) {
  mean(intensity)
}

findChromPeaks <- function(spec, sn=5, rtgap=0.2 ) {

  noise <- estimateChromPeakNoise(spec[,"intensity"] )

  spectab <- matrix(nrow = 0, ncol = 4)
  colnames(spectab) <- c("rt", "rt.min", "rt.max", "sn")

  while (spec[i <- which.max(spec[,"intensity"]), "intensity"] > noise*sn) {

    rt <- spec[i,"rt"]
    intensity <- spec[i,"intensity"]
    rt.range <- xcms:::descendValue(spec[,"intensity"], noise, i)
    rt.range <- c(max(1,rt.range[1]-1), min(nrow(spec),rt.range[2]+1) )

    if (rt.range[1] >= 1 && rt.range[2] <= nrow(spec)) {
      rt.min <- spec[rt.range[1],"rt"]
      rt.max <- spec[rt.range[2],"rt"]
      if (!any(abs(spectab[,"rt"] - rt) <= rtgap))
        spectab <- rbind(spectab, c(rt, rt.min, rt.max, spec[i,"intensity"]/noise))
    }
    spec[seq(rt.range[1], rt.range[2]),"intensity"] <- 0
    }

    spectab
}

findPairChromPeaks <- function(rt, light.int, heavy.int, rt.range, sn=5) {
  rtdiff <- mean( diff(rt) )

  m.light <- cbind(rt,light.int)
  dimnames(m.light)[[2]] <- c("rt","intensity")
  m.light <- m.light[rt>=rt.range[1]&rt<=rt.range[2],]
  noise.light <- estimateChromPeakNoise(m.light[,"intensity"])
  peaks.light <- findChromPeaks(m.light,sn,rtgap=0.5)
  n.light <- dim(peaks.light)[[1]]

  m.heavy <- cbind(rt,heavy.int)
  dimnames(m.heavy)[[2]] <- c("rt","intensity")
  m.heavy <- m.heavy[rt>=rt.range[1]&rt<=rt.range[2],]
  noise.heavy <- estimateChromPeakNoise(m.heavy[,"intensity"])
  peaks.heavy <- findChromPeaks(m.heavy,sn,rtgap=5*rtdiff)
  n.heavy <- dim(peaks.heavy)[[1]]

  pair.range <- c(noise.light, noise.heavy)
  if ( n.heavy == 0 | n.light == 0 ) return(pair.range)

  for (i in 1:n.light) {
    rt.i <- peaks.light[i,"rt"]
    d.rt <- abs(rt.i-peaks.heavy[,"rt"])
    j <- which.min(d.rt)
   # if (d.rt[j]>5*rtdiff) next
    rt.j <- peaks.heavy[j,"rt"]
    if ( rt.j >=peaks.light[i,"rt.min"] & rt.j <= peaks.light[i,"rt.max"]
        &  rt.i >=peaks.heavy[j,"rt.min"] & rt.i <= peaks.heavy[j,"rt.max"] ) {
      low <- min(peaks.light[i,"rt.min"],peaks.heavy[j,"rt.min"])
      high <- max(peaks.light[i,"rt.max"],peaks.heavy[j,"rt.max"])
      if ( low < high ) {
        pair.range <- c(pair.range,low,high)
      }
    }
  }
  return(pair.range)
}

#findPairChromPeaks <- function(rt, light.int, heavy.int, rt.range) {
#  m.light <- cbind(rt,light.int)
#  dimnames(m.light)[[2]] <- c("mz","intensity")
#  m.light <- m.light[rt>=rt.range[1]&rt<=rt.range[2],]
#  noise.light <- specNoise(m.light)
#  peaks.light <- specPeaks(m.light,sn=10,mzgap=0.5)
#  n.light <- dim(peaks.light)[[1]]

#  m.heavy <- cbind(rt,heavy.int)
#  dimnames(m.heavy)[[2]] <- c("mz","intensity")
#  m.heavy <- m.heavy[rt>=rt.range[1]&rt<=rt.range[2],]
#  noise.heavy <- specNoise(m.heavy)
#  peaks.heavy <- specPeaks(m.heavy,sn=10,mzgap=0.5)
#  n.heavy <- dim(peaks.heavy)[[1]]

#  pair.range <- integer(0)
#  if ( n.heavy == 0 | n.light == 0 ) return(pair.range)

#  for (i in 1:n.light) {
#    rt.i <- peaks.light[i,"mz"]
#    d.rt <- abs(rt.i-peaks.heavy[,"mz"])
#    j <- which.min(d.rt)
#    if (d.rt[j]>0.5) next
#    rt.j <- peaks.heavy[j,"mz"]
#    low <- min(rt.i-1.2*peaks.light[i,"fwhm"],rt.j-1.2*peaks.heavy[j,"fwhm"])
#    high <- max(rt.i+1.2*peaks.light[i,"fwhm"],rt.j+1.2*peaks.heavy[j,"fwhm"])
#    pair.range <- c(pair.range,low,high)
#  }
#  return(pair.range)
#}

readFileFromMsn <- function( xcms.raw ) {
  filename <- xcms.raw@filepath
  filename.base <- strsplit(filename,".mzXML")[[1]][1]
  dtable <- read.table(paste(filename.base,".ms1.mz_int",sep=""), header=T)
  xcms.raw@env$mz <- dtable[,"mz"]
  xcms.raw@env$intensity <- dtable[,"int"]
  dtable <- read.table(paste(filename.base,".ms1.index_table",sep=""), header=T)
  xcms.raw@scanindex <- dtable[,"scanindex"]
  xcms.raw@scantime  <- dtable[,"scantime"]
  xcms.raw@acquisitionNum <- dtable[,"acquisitionNum"]
  if ( length(xcms.raw@msnScanindex) > 0 ) {
    dtable <- read.table(paste(filename.base,".ms2.mz_int",sep=""), header=T)
    xcms.raw@env$msnMz <- dtable[,"mz"]
    xcms.raw@env$msnIntensity <- dtable[,"int"]
    dtable <- read.table(paste(filename.base,".ms2.index_table",sep=""), header=T)
    xcms.raw@msnScanindex <- dtable[, "msnScanindex"]
    xcms.raw@msnAcquisitionNum <- dtable[, "msnAcquisitionNum"]
    xcms.raw@msnPrecursorScan <- dtable[, "msnPrecursorScan"]
    xcms.raw@msnLevel <- dtable[, "msnLevel"]
    xcms.raw@msnRt <- dtable[, "msnRt"]
    xcms.raw@msnPrecursorMz <- dtable[, "msnPrecursorMz"]
    xcms.raw@msnPrecursorIntensity <- dtable[, "msnPrecursorIntensity"]
    xcms.raw@msnPrecursorCharge <- rep(1, length(xcms.raw@msnPrecursorIntensity))
    xcms.raw@msnCollisionEnergy <- rep(xcms.raw@msnCollisionEnergy[1], length(xcms.raw@msnPrecursorIntensity))
  }
  return( xcms.raw )
}
