library(xcms)

#memory.limit(2560) # 2.5G memory maximum

input.file <- "bash.input.file"

input.file.suffix <- ".mzXML"
output.file.suffix <- ".png"

# within 10sec of retention window
rt.cut <- 10.0
# ppm
mz.ppm.cut <- 0.000025 # 25ppm
# From Eranthie's isotopically labeled probe
pair.mass.delta <- 6.01381
# nature mass difference between C12 and C13
isotope.mass.unit <- 1.0033548
mass.shift <- round( pair.mass.delta/isotope.mass.unit )
# mass of a proton
Hplus.mass <- 1.0072765
# possible charge states
charge.states <- seq(1,6)
# peak intensity ratio
intensity.tag <- "maxo"
# ratio of intensity of heavy and light peaks
intensity.ratio.range <- c(0.5,2.0)

# find all matched input files in current directory
input.path <- getwd()
#input.files <- list.files( path=input.path, pattern=paste(input.file.suffix,"$",sep="") )

#for ( input.file in input.files ) {
input.file.base <- strsplit( input.file, input.file.suffix )[[1]][1]
output.file.base <- paste(input.file.base,"doublet",sep="_")
output.path <- paste(input.path,"/", input.file.base, "_doublet_peaks/",sep="")
dir.create(output.path)

filename <- paste(input.path, "/", input.file.base, input.file.suffix, sep="")
xfile <- xcmsRaw( filename, profstep=0 )

xpeaks <- findPeaks(xfile,method="centWave", ppm=5, snthresh=2.0, minptsperpeak=2)
                                        # sort the matrix by retention time first
xpeaks.save <- xpeaks
xpeaks <- xpeaks[order(xpeaks[,"rt"]), ]

                                        # label rownames
num.peaks <- dim(xpeaks)[1]
rownames(xpeaks) <- seq(1, num.peaks)

# table with raw pair peak index, charge, rt and hit.id etc
out.filename0 <- paste(output.path, output.file.base, ".xcms_peaks.txt", sep="")
write.table( format(xpeaks,digits=10), out.filename0, quote=FALSE, row.names=FALSE)

                                        # make a fresh empty list
n.nonredundant.hit <- 0
pair.list <- data.frame(idx1=numeric(0),idx2=numeric(0),charge=numeric(0),
                        mz1=numeric(0), mz2=numeric(0),
                        rt1=numeric(0), rt2=numeric(0), hit.id=numeric(0) )

rt.list <- levels(factor(xpeaks[,"rt"]))
cat( paste("detect isotopic peak pairs in", length(rt.list), "rt windows in total\n") )
cat( "Finished " )
for ( rt.i in 1:length(rt.list) ) {
  cat(".")
  if (rt.i%%100 == 0) cat(rt.i)
  rt <- rt.list[rt.i]
  rt.range <- c( as.numeric(rt)-rt.cut, as.numeric(rt)+rt.cut)
  this.rt <- (xpeaks[,"rt"] > rt.range[1]) & (xpeaks[,"rt"] < rt.range[2] )
  #this.rt <- factor(xpeaks[,"rt"])==rt
  if ( sum(this.rt) <= 1 ) next
                                        # get all peaks with this retension time and sort by mz
  xpeaks.this.rt <- xpeaks[this.rt,]
  xpeaks.this.rt <- xpeaks.this.rt[order(xpeaks.this.rt[,"mz"]), ]
  mz.this.rt <- xpeaks.this.rt[,"mz"]
  num.mz.this.rt <- length( mz.this.rt )
  mz.delta <- matrix(0, num.mz.this.rt, num.mz.this.rt,
                     dimnames=list(names(mz.this.rt), names(mz.this.rt) ) )
  mz.delta.ppm <- mz.delta # tolerance matrix

  # for redundancy check
  # pair.list.prev.rt contains all existing pairs
  # which could potentially overlap with the current rt window
  pair.list.prev.rt <- pair.list[ (pair.list[,"rt1"] >= rt.range[1] | pair.list[,"rt2"] >= rt.range[1] ), ]

  for (i in 1:num.mz.this.rt ) {
    for ( j in i:num.mz.this.rt ) {
      mz.delta[i,j] <- mz.delta[j,i] <- abs(mz.this.rt[i] - mz.this.rt[j])
      mz.delta.ppm[i,j] <- mz.delta.ppm[j,i]  <- mz.ppm.cut * ( mz.this.rt[i]+mz.this.rt[j] ) / 2.0
    }
  }
  for (charge in charge.states) {
    pair.list.this.charge <- pair.list[F,] # an empty dataframe tracking pairs within this rt and charge
    n.hit.this.charge <- 0
    pair.mz.delta <- pair.mass.delta / charge
    mz.delta.pass <- abs(mz.delta - pair.mz.delta) < mz.delta.ppm
    for (i in 1:num.mz.this.rt ) {
      for ( j in i:num.mz.this.rt ) {
        if ( mz.delta.pass[i,j] ) {
                                        # earlier sorting by mz ensures that 1 is alway lighter peak
          peak.index1 <- dimnames(mz.delta.pass)[[1]][i]
          peak.index2 <- dimnames(mz.delta.pass)[[2]][j]
                                        # check proper isotopic distribution given this charge
                                        # simple criteria -- for each of light and heavy peaks,
                                        # there is at least one isotopic peak nearby
          isotope.mz.delta <- isotope.mass.unit / charge
          isotope.delta.pass <- abs(mz.delta - isotope.mz.delta) < mz.delta.ppm
          if ( (sum(isotope.delta.pass[i,]) == 0) | (sum(isotope.delta.pass[,j]) == 0) ) next
                                        # mz
          peak.mz1 <- xpeaks[peak.index1, "mz"]
          peak.mz2 <- xpeaks[peak.index2, "mz"]
                                        # intensity
          peak.maxo1 <- xpeaks[peak.index1, intensity.tag]
          peak.maxo2 <- xpeaks[peak.index2, intensity.tag]
                                        # retension time
          peak.rt1 <- xpeaks[peak.index1, "rt"]
          peak.rt2 <- xpeaks[peak.index2, "rt"]
                                        # check peak intensity ratio to fall in certain range
          ratio <- peak.maxo1 / peak.maxo2
          if ( ratio < intensity.ratio.range[1] | ratio > intensity.ratio.range[2] ) next

          # check whether the exact peak pair has been found in previous rt window before
          redundant <- FALSE
          if ( dim(pair.list.prev.rt)[1] > 0 ) {
            for (k in 1:dim(pair.list.prev.rt)[1] ) {
              if ((pair.list.prev.rt[k,"idx1"] == peak.index1) &&
                  (pair.list.prev.rt[k,"idx2"] == peak.index2) &&
                  (pair.list.prev.rt[k,"charge"] == charge ) ) {
                redundant <- TRUE
                break
              }
            }
          }
          if ( redundant ) next

          # this pair is uniq, but maybe it belongs an existing isotope envelope pair
          # the last digit is used to track unique isotope envelope pairs
          new.pair <- pair.list.this.charge[F,]
          #new.pair[1,] <- c( peak.index1, peak.index2, charge, max(peak.rt1, peak.rt2), 0 )
          new.pair[1,] <- c( peak.index1, peak.index2, charge, peak.mz1, peak.mz2, peak.rt1, peak.rt2, 0 )
          nrow.this.charge <- dim(pair.list.this.charge)[1]
          if ( nrow.this.charge == 0 ) {
            n.hit.this.charge <- n.hit.this.charge + 1
            new.pair[1,"hit.id"] <- n.hit.this.charge
            pair.list.this.charge <- rbind( pair.list.this.charge, new.pair )
          } else {
            for ( k in 1:nrow.this.charge) {
              tmp.mz1 <- xpeaks[ pair.list.this.charge[k,"idx1"], "mz" ]
              tmp.mz2 <- xpeaks[ pair.list.this.charge[k,"idx2"], "mz" ]
              if (((tmp.mz2>peak.mz1)&(tmp.mz2<peak.mz2)) |
                  ((tmp.mz1>peak.mz1)&(tmp.mz1<peak.mz2)) ) {
                # this pair is from an existing isotopic envelope,
                new.pair[1,"hit.id"] <- pair.list.this.charge[k,"hit.id"]
                break
              }
            }
            # this pair belongs to a new pair of isotopic envelopes
            if ( new.pair[1,"hit.id"] == 0 ) {
              n.hit.this.charge <- n.hit.this.charge + 1
              new.pair[1,"hit.id"] = n.hit.this.charge
            }
            # add this pair to the list
            pair.list.this.charge <- rbind( pair.list.this.charge, new.pair )
          }
        }
      }
    }
    if (n.hit.this.charge == 0) next
    pair.list.tmp <- pair.list.this.charge[F,]
    for( k in 1:n.hit.this.charge ) {
      pair.list.tmp <- pair.list.this.charge[pair.list.this.charge[,"hit.id"] == k, ]
      rt.factor <- factor( c(pair.list.tmp$rt1,pair.list.tmp$rt2) )
      rt.levels <- levels(rt.factor)
      main.rt <- rt.levels[ order(table(rt.factor),rev(rt.levels))[length(rt.levels)] ]
      if ( main.rt == rt ) {
        n.nonredundant.hit <- n.nonredundant.hit + 1
        pair.list.tmp[,"hit.id"] <- n.nonredundant.hit
        pair.list<- rbind( pair.list, pair.list.tmp )
      }
    }
    #pair.list <- rbind( pair.list, pair.list.this.charge )
  }
}
cat("\n")

if( dim(pair.list)[1] == 0 ) {
  print("no doublet peaks identified")
  quit()
}
# table with raw pair peak index, charge, rt and hit.id etc
out.filename0 <- paste(output.path, output.file.base, ".pair_index.txt", sep="")
write.table( format(pair.list,digits=10), out.filename0, quote=FALSE, row.names=FALSE)

# filter and output
rt.window.per.scan <- 6.0 # time elapse per ms1 scan
mass.window <- 2.0  # plot 2.0 m/z unit from center of delected peak
num.ms2.per.ms1 <- 18
out.table <- data.frame()
for ( id in levels(factor(pair.list[,"hit.id"])) ) {
  pair.list.this.id <- pair.list[ pair.list[,"hit.id"]==id, ]
  charge <- as.numeric(pair.list.this.id[1,"charge"])
  isotope.mz.delta <- isotope.mass.unit / charge
  ## peak intensity
  peaks.maxo1 <- xpeaks[ as.character( pair.list.this.id[,"idx1"] ), "maxo"]
  peaks.maxo2 <- xpeaks[ as.character( pair.list.this.id[,"idx2"] ), "maxo"]
  peaks.maxo <- peaks.maxo1 + peaks.maxo2
  peaks.maxo.cut <- peaks.maxo > max(peaks.maxo) * 0.1
  ## within same isotopic cluster, get rid background noise by 1% cut
  peaks.maxo <- peaks.maxo[peaks.maxo.cut]
  peaks.maxo1 <- peaks.maxo1[peaks.maxo.cut]
  peaks.maxo2 <- peaks.maxo2[peaks.maxo.cut]
  pair.list.this.id <- pair.list.this.id[peaks.maxo.cut,]
  peaks.mz1 <- as.numeric(pair.list.this.id[,"mz1"])
  peaks.mz2 <- as.numeric(pair.list.this.id[,"mz2"])
  ## most intense peak position in the light envelope
  i <- which.max(peaks.maxo1)
  envelope1 <- abs( (peaks.mz1-peaks.mz1[i])/isotope.mz.delta )
  envelope2 <- abs( (peaks.mz2-peaks.mz2[i])/isotope.mz.delta )
  envelope.check <- (abs(envelope1-round(envelope1))/peaks.mz1[i] < mz.ppm.cut &
                     abs(envelope2-round(envelope2))/peaks.mz2[i] < mz.ppm.cut &
                     round(envelope1) <= 4 &
                     round(envelope2) <= 4)
  pair.list.this.id <- pair.list.this.id[envelope.check,]
  peaks.mz1 <- as.numeric(pair.list.this.id[,"mz1"])
  peaks.mz2 <- as.numeric(pair.list.this.id[,"mz2"])
  peaks.maxo <- peaks.maxo[envelope.check]
  peaks.maxo1 <- peaks.maxo1[envelope.check]
  peaks.maxo2 <- peaks.maxo2[envelope.check]

  redundancy.check <- which( diff(peaks.mz1)/peaks.mz1[-1] < mz.ppm.cut ) + 1
  if ( length(redundancy.check)>0 ) {
    pair.list.this.id <- pair.list.this.id[-redundancy.check,]
    peaks.mz1 <- as.numeric(pair.list.this.id[,"mz1"])
    peaks.mz2 <- as.numeric(pair.list.this.id[,"mz2"])
    peaks.maxo <- peaks.maxo[-redundancy.check]
    peaks.maxo1 <- peaks.maxo1[-redundancy.check]
    peaks.maxo2 <- peaks.maxo2[-redundancy.check]
  }

  ## num of peaks
  peaks.n <- length( peaks.maxo )
  ## need at least three peak pairs per hit
  if (peaks.n < 3) next

  ## calculate isotope distribution based on averagine assumption
  predicted.mass <- charge*peaks.mz1[1]
  predicted.comp <- averagine.count( predicted.mass)
  predicted.dist <- isotope.dist(predicted.comp)
  predicted.dist.merge <- c(predicted.dist,rep(0,mass.shift)) + c(rep(0,mass.shift),predicted.dist)
  predicted.dist.merge <- predicted.dist.merge/max(predicted.dist.merge)
  observed.mz <- c(peaks.mz1, peaks.mz2)
  observed.maxo <- c(peaks.maxo1,peaks.maxo2)
  observed.maxo <- observed.maxo/max(observed.maxo)
  observed.index <- round((observed.mz - peaks.mz1[1])/isotope.mz.delta)
  pcor <- rep(0,mass.shift)
  for (s in 1:mass.shift) {
    predicted.maxo <- predicted.dist.merge[observed.index+s]
    pcor[s] <-  lsfit(predicted.maxo,observed.maxo)[[1]][2]#cov(predicted.maxo,observed.maxo)/(sd(predicted.maxo)*sd(observed.maxo))
  }
  ##best.shift <- which.max(pcor) - 1
  ##best.pcor <- max(pcor)
  best.shift <- which.min(abs(pcor-1.0))
  best.pcor <- pcor[best.shift]
  best.shift <- best.shift-1
  predicted.mono.mz <- peaks.mz1[1] - best.shift*isotope.mz.delta
  upper.bound <- max( which(predicted.dist.merge> 0.01) )
  upper.bound.mz <- predicted.mono.mz + (upper.bound-1)*isotope.mz.delta
  predicted.mz <- predicted.mono.mz + seq(0,upper.bound-1)*isotope.mz.delta
  predicted.int <- predicted.dist.merge[1:upper.bound]
  ## plot chromatogram of the most intense peak,
  ## label others under the same envelope
  ## most intense peak position in the light envelope
  i <- which.max(peaks.maxo1)
  rt <- as.numeric( pair.list.this.id[i,"rt1"] )
  ##rt <- ( as.numeric(pair.list.this.id[i,"rt1"]) + as.numeric(pair.list.this.id[i,"rt2"]) ) / 2.0
  rt.range <- c( rt-rt.window.per.scan, rt+rt.window.per.scan)
  ## get corresponding scan number
  best.scan.number <- which(xfile@scantime == rt)
  ## for plot a specific ms1 scan
  mass.range <- c( min(peaks.mz1) - mass.window, max(peaks.mz2) + mass.window )
  ## get actual ms1 scan spectrum
  scan.data <- getScan(xfile, best.scan.number, massrange=mass.range)
  scan.mz <- scan.data[,"mz"]
  scan.int <- scan.data[,"intensity"]
  tmp <- (scan.mz >= min(peaks.mz1) & scan.mz <= max(peaks.mz2) )
  ##scan.mz.zoom <- scan.mz[ tmp ]
  scan.int.zoom <- scan.int[ tmp ]
  ## the most intensive peak should match one of the pairs
  ## tmp.mz <- scan.mz.zoom[ order(scan.int.zoom)[length(scan.int.zoom)] ]
  if ( max(peaks.maxo1,peaks.maxo2) < 0.9*max(scan.int.zoom) ) next
  ## the most intensive peak
  peak.mz1 <- peaks.mz1[i]
  peak.mz2 <- peaks.mz2[i]
  peak.maxo1 <- scan.int[ abs(scan.mz-peak.mz1) < mz.ppm.cut * peak.mz1 ]
  peak.maxo2 <- scan.int[ abs(scan.mz-peak.mz2) < mz.ppm.cut * peak.mz2 ]
  ##overlapping peaks
  if ( length(peak.maxo1) != 1 | length(peak.maxo2) != 1 ) next

  out.filename <- paste(output.path, output.file.base, "_", peaks.n, "_", id, output.file.suffix, sep="")
  ## make plot for each isotopic cluster
  png(out.filename)
  par(mfrow=c(2,1))
  par(las=0)
  ## mz spectrum top
  ##plotScan(xfile, best.scan.number, massrange=mass.range)
  scan.int <- scan.int/max(scan.int.zoom)
  ##least-sqaured fit with single scan ###
  observed.int <- rep(0,upper.bound)
  for ( i in 1:upper.bound ) {
    i.mz <- predicted.mz[i]
    matched.index <- which(abs(scan.mz-i.mz)/predicted.mono.mz<mz.ppm.cut)
    if (length(matched.index)>1) {
      matched.index <- matched.index[ which.min(abs(scan.mz[matched.index]-predicted.mono.mz)) ]
    }
    i.int <- scan.int[matched.index]
    if (length(i.int) == 0) {
      i.int = 0
    }
    observed.int[i] <- i.int
  }
  observed.mono.int <- observed.int[1]
  best.pcor <-  lsfit(predicted.int,observed.int)[[1]][2] #coefficient-slope
  ## plot ###
  ylimit <- c( 0, 1.05) #max(scan.int[(scan.mz>=predicted.mono.mz&scan.mz<=upper.bound.mz)]) )
  plot( scan.mz, scan.int, type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit)

  title( paste("Scan # ", best.scan.number, " @ ", round(xfile@scantime[best.scan.number], 1),
               " sec (raw # ", (best.scan.number-1)*(num.ms2.per.ms1+1)+1, " @ ",
               round(xfile@scantime[best.scan.number]/60,1)," min)", sep = "")
        )
  title( paste("Charge:",charge, "; Mono.mz:", round(predicted.mono.mz,3),"(M); NL:",
               formatC(max(peak.maxo1,peak.maxo2), digits=2,format="e")), line=0.5)
  ##ymark <- min( scan.data[,2] ) + 10
  ymark <- 0
  ##label all identified peaks under this isotopic envelope
  points( peaks.mz1, rep(ymark, length(peaks.mz1) ), pch=23, col="red",bg="red")
  points( peaks.mz2, rep(ymark, length(peaks.mz2) ), pch=24, col="blue",bg="blue")
  ##points( peak.mz1, peak.maxo1, pch="C", col="red",bg="red")
  ##points( peak.mz2, peak.maxo2, pch="C", col="blue",bg="blue")
  points( predicted.mono.mz, min(observed.mono.int,1.0), pch="M", col="red",bg="red")
  par(new=T)
  plot( predicted.mz, predicted.int, type='b',xlab="",ylab="",col="green",axes=F,xlim=mass.range,ylim=ylimit)
  ##chromatogram bottom
  raw.ECI.light <- rawEIC(xfile, c(peak.mz1*(1-mz.ppm.cut), peak.mz1*(1+mz.ppm.cut)) )
  raw.ECI.heavy <- rawEIC(xfile, c(peak.mz2*(1-mz.ppm.cut), peak.mz2*(1+mz.ppm.cut)) )
  xlimit <-c( max(1, best.scan.number-5), min(best.scan.number+5, length(raw.ECI.light[[1]]) ) )
  ylimit <- range(c(raw.ECI.light[[2]][xlimit[1]:xlimit[2]], raw.ECI.heavy[[2]][xlimit[1]:xlimit[2]]))
  plot(raw.ECI.light[[1]], raw.ECI.light[[2]], type="l", col="red",xlab="scan #", ylab="intensity",
       main=paste("Chromatogram of", as.character(format(peak.mz1, digits=7)),
         "and", as.character(format(peak.mz2,digits=7)), "m/z (C)"), xlim=xlimit, ylim=ylimit)
  lines(raw.ECI.heavy[[1]], raw.ECI.heavy[[2]], col='blue', xlim=xlimit, ylim=ylimit)
  dev.off()
  ## save entries in a text table for output
  outlist <- list( best.scan.number, charge, predicted.mono.mz*charge-Hplus.mass*charge, predicted.mono.mz, peak.mz1, as.numeric(id), peaks.n, best.pcor )
  names(outlist) <- c( "scan", "charge", "mono.mass", "mono.mz", "peak.mz", "peaks.id", "peaks.count", "pcor")
  if ( length(out.table) == 0 ) {
    out.table <- data.frame( outlist )
  } else {
    out.table <- rbind( out.table, outlist )
  }
}

## table with peak.count
rownames(out.table) <- out.table$peaks.id
out.table <- out.table[ order(out.table$peaks.count, decreasing=TRUE), ]
out.filename <- paste(output.path, output.file.base, ".table.txt", sep="")
write.table( format(out.table,digits=10), out.filename, quote=FALSE, row.names=FALSE)

## make a pep3d-style xc/ms intensity plot
out.filename2 <- paste(output.path, output.file.base, "_MS_SCAN", output.file.suffix, sep="")
png(out.filename2,width=8.5,height=11,units="in",res=72)
par(mfrow=c(2,1))
par(las=1)

prof.range <- profRange(xfile)
masslim <- round( prof.range$massrange )
scanlim <- prof.range$scanrange
massname <- masslim[1] : masslim[2]
scanname <- scanlim[1] : scanlim[2]
prof.matrix <- matrix( 0.0,nrow=length(massname), ncol=length(scanname), dimnames=list(massname,scanname) )

for (i in 1:length(xfile@scanindex)) {
  idx <- (xfile@scanindex[i]+1):min(xfile@scanindex[i+1],length(xfile@env$mz), na.rm=TRUE)
  prof.matrix[as.character(round(xfile@env$mz[idx])), as.character(i)] <- xfile@env$intensity[idx]
}
image(massname, xfile@scantime[scanname], log(prof.matrix+1), col=gray(256:0/256),
      xlab="m/z", ylab="Seconds", zlim=log(range(xfile@env$intensity)) )
title("XC/MS Log Intensity Image (whole)")
col.max <- max( out.table[,"peaks.count"] )
for (i in 1:dim(out.table)[1] ) {
  xp <- round( out.table[i,"mono.mz"] )
  yp <- out.table[i,"scan"]
  countp <- out.table[i,"peaks.count"]
  points( xp, xfile@scantime[yp], type="o", col = rainbow(col.max)[countp] )
}

x.lim <- round (range(out.table[,"mono.mz"]) )
y.lim <- range(xfile@scantime[out.table[,"scan"]])
image(massname, xfile@scantime[scanname], log(prof.matrix+1), col=gray(256:0/256),
      xlab="m/z", ylab="Seconds", zlim=log(range(xfile@env$intensity)), xlim=x.lim, ylim=y.lim )
title("XC/MS Log Intensity Image (zoom)")
col.max <- max( out.table[,"peaks.count"] )
for (i in 1:dim(out.table)[1] ) {
  xp <- round( out.table[i,"mono.mz"] )
  yp <- out.table[i,"scan"]
  countp <- out.table[i,"peaks.count"]
  points( xp, xfile@scantime[yp], type="o", col = rainbow(col.max)[countp] )
}
dev.off()

