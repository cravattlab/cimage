library(xcms)

input.path <- "C:/cygwin/home/chuwang/light/"
input.file.base <- "test_cen"
input.file.suffix <-".mzXML"
output.path <- paste(input.path,"xcms_cen_plot/",sep="")
output.file.base <- paste(input.file.base,"doublet",sep="_")
output.file.suffix <- ".png"

filename <- paste(input.path, input.file.base, input.file.suffix, sep="")
xfile <- xcmsRaw( filename, profstep=0 )

xpeaks <- findPeaks(xfile,method="centWave", ppm=25, snthresh=2.0,minptsperpeak=5)
num.peaks <- dim(xpeaks)[1]
rownames(xpeaks) <- seq(1, num.peaks)

# sort the matrix by retention time first
xpeaks <- xpeaks[order(xpeaks[,"rt"]), ]

# within 10sec of retention window
rt.cut <- 10.0

mz.ppm.cut <- 0.000025 # 25ppm
# From Eranthie's isotopically labeled probe
pair.mass.delta <- 6.01381
# nature mass difference between C12 and C13
isotope.mass.unit <- 1.00286864
# possible charge states
charge.states <- seq(1,5)
# peak intensity ratio
intensity.tag <- "maxo" #
intensity.ratio.range <- c(0.5,2.0) # ratio of intensity of heavy and light peaks

# make a fresh empty list
n.nonredundant.hit <- 0
pair.list <- NULL

for ( rt in levels(factor(xpeaks[,"rt"])) ) {
#  rt.range <- c( as.numeric(rt)-rt.cut, as.numeric(rt)+rt.cut)
#  this.rt <- (xpeaks[,"rt"] > rt.range[1]) & (xpeaks[,"rt"] < rt.range[2] )
  this.rt <- factor(xpeaks[,"rt"])==rt
  if ( sum(this.rt) <= 1 ) next
  # get all peaks with this retension time and sort by mz
  xpeaks.this.rt <- xpeaks[this.rt,]
  xpeaks.this.rt <- xpeaks.this.rt[order(xpeaks.this.rt[,"mz"]), ]
  mz.this.rt <- xpeaks.this.rt[,"mz"]
  num.mz.this.rt <- length( mz.this.rt )
  mz.delta <- matrix(0, num.mz.this.rt, num.mz.this.rt, dimnames=list(names(mz.this.rt), names(mz.this.rt) ) )
  mz.delta.ppm <- mz.delta # tolerance matrix

  for (i in 1:num.mz.this.rt ) {
    for ( j in i:num.mz.this.rt ) {
      mz.delta[i,j] <- mz.delta[j,i] <- abs(mz.this.rt[i] - mz.this.rt[j])
      mz.delta.ppm[i,j] <- mz.delta.ppm[j,i]  <- mz.ppm.cut * ( mz.this.rt[i]+mz.this.rt[j] ) / 2.0
    }
  }
  for (charge in charge.states) {
    pair.list.this.charge <- NULL
    pair.mz.delta <- pair.mass.delta / charge
    mz.delta.pass <- abs(mz.delta - pair.mz.delta) < mz.delta.ppm
    for (i in 1:num.mz.this.rt ) {
      for ( j in i:num.mz.this.rt ) {
        if ( mz.delta.pass[i,j] ) {
          # earlier sorting by mz ensures that 1 is alway lighter peak
          peak.index1 <- dimnames(mz.delta.pass)[[1]][i]
          peak.index2 <- dimnames(mz.delta.pass)[[2]][j]
          # check proper isotopic distribution given this charge
          # simple criteria -- for each of light and heavy peaks, there is at least one isotopic peak nearby
          isotope.mz.delta <- isotope.mass.unit / charge
          isotope.delta.pass <- abs(mz.delta - isotope.mz.delta) < mz.delta.ppm
          if ( (sum(isotope.delta.pass[i,]) == 0) | (sum(isotope.delta.pass[,j]) == 0) ) next
          # mz
          peak.mz1 <- xpeaks[peak.index1, "mz"]
          peak.mz2 <- xpeaks[peak.index2, "mz"]
          # intensity
          peak.maxo1 <- xpeaks[peak.index1, intensity.tag]
          peak.maxo2 <- xpeaks[peak.index2, intensity.tag]
          #check peak intensity ratio to fall in certain range
          ratio <- peak.maxo1 / peak.maxo2
          if ( ratio < intensity.ratio.range[1] | ratio > intensity.ratio.range[2] ) next
          new.pair <- c( peak.index1, peak.index2, charge, "0" ) # last number is tracking nonredudant hit
          if ( length(pair.list.this.charge) == 0 ) {
            n.nonredundant.hit <- n.nonredundant.hit + 1
            new.pair[4] <- n.nonredundant.hit
            pair.list.this.charge <- new.pair
          } else {
            for ( k in seq(1,length(pair.list.this.charge),by=4) ) {
              tmp.mz1 <- xpeaks[ pair.list.this.charge[k], "mz" ]
              tmp.mz2 <- xpeaks[ pair.list.this.charge[k+1], "mz" ]
              if ( ((tmp.mz2>peak.mz1)&(tmp.mz2<peak.mz2)) | ((tmp.mz1>peak.mz1)&(tmp.mz1<peak.mz2)) ) {
                # this pair is from an existing isotopic cluster
                new.pair[4] <- pair.list.this.charge[k+3]
                break
              }
            }
            if ( new.pair[4] == "0" ) {
              n.nonredundant.hit <- n.nonredundant.hit + 1
              new.pair[4] = n.nonredundant.hit
            }
            pair.list.this.charge <- c( pair.list.this.charge, new.pair)
          }
        }
      }
    }
    if ( length(pair.list) == 0 ) {
      pair.list <- pair.list.this.charge
    } else {
      pair.list <- c( pair.list, pair.list.this.charge )
    }
    pair.list.this.charge <- NULL
  }
}
pair.matrix <- matrix(pair.list, ncol=4, byrow=TRUE)
colnames(pair.matrix) <- c("idx1", "idx2", "charge", "hit.id")
pair.matrix <- data.frame( pair.matrix )

rt.window.per.scan <- 6.0 # time elapse per ms1 scan
mass.window <- 2.0 # plot 2.0 m/z unit from center of delected peak
num.ms2.per.ms1 <- 18
for ( id in levels(factor(pair.matrix[,"hit.id"])) ) {
  pair.matrix.this.id <- pair.matrix[ pair.matrix[,"hit.id"]==id, ]
  if ( dim(pair.matrix.this.id)[1] > 1) {
    prefix <- "M_"
  } else {
    prefix <- "S_"
  }
  peaks.mz1 <- xpeaks[ as.character( pair.matrix.this.id[,"idx1"] ), "mz"]
  peaks.mz2 <- xpeaks[ as.character( pair.matrix.this.id[,"idx2"] ), "mz"]
  # for plot a specific ms1 scan
  mass.range <- c( min(peaks.mz1) - mass.window, max(peaks.mz2) + mass.window )
  # for peak intensity
  peaks.maxo1 <- xpeaks[ as.character( pair.matrix.this.id[,"idx1"] ), "maxo"]
  peaks.maxo2 <- xpeaks[ as.character( pair.matrix.this.id[,"idx2"] ), "maxo"]
  peaks.maxo <- peaks.maxo1 + peaks.maxo2
  # num of peaks
  peaks.n <- length( peaks.maxo )
  # plot chromatogram of the most intense peak, label others under the same envelope
  i <- order(peaks.maxo)[length(peaks.maxo)]
  peak.index1 <- as.character(pair.matrix.this.id[i,1])
  peak.index2 <- as.character(pair.matrix.this.id[i,2])
  charge <- pair.matrix.this.id[i,3]

  rt <- ( xpeaks[peak.index1,"rt"] + xpeaks[peak.index2,"rt"] ) / 2.0
  rt.range <- c( rt-rt.window.per.scan, rt+rt.window.per.scan)

  raw.EIC.data <- rawEIC(xfile, massrange=mass.range, timerange=rt.range)
  best.scan.number <- raw.EIC.data$scan[ order(raw.EIC.data$intensity)[length(raw.EIC.data$intensity)] ]
  # get actual ms1 scan spectrum
  scan.data <- getScan(xfile, best.scan.number, massrange=mass.range)
  scan.mz <- scan.data[,"mz"]
  scan.int <- scan.data[,"intensity"]
  peak.mz1 <- peaks.mz1[i]
  peak.mz2 <- peaks.mz2[i]
  peak.maxo1 <- scan.int[ abs(scan.mz-peak.mz1) < mz.ppm.cut * peak.mz1 ]
  peak.maxo2 <- scan.int[ abs(scan.mz-peak.mz2) < mz.ppm.cut * peak.mz2 ]
  # a very simple way to find the monoisotopic peak
  isotope.mz.delta <- isotope.mass.unit / as.numeric(as.character(charge))
  lower.mz <- monoiso.mz <- peak.mz1
  lower.int <- monoiso.int <- peak.maxo1
  lower.int.ratio <- 1.0
  repeat {
    lower.mz.exist <- abs(scan.mz - monoiso.mz + isotope.mz.delta)/monoiso.mz < mz.ppm.cut
    if ( sum( lower.mz.exist ) < 1 | lower.int.ratio < 0.5) {
      # cannot find lower peak any more or intensity not strong enough, finish searching
      break
    } else {
      monoiso.mz <- lower.mz
      monoiso.int <- lower.int
      lower.mz <- scan.mz[lower.mz.exist][1]
      lower.int <- scan.int[lower.mz.exist][1]
      lower.int.ratio <- lower.int/monoiso.int
    }
  } # repeat
  # simililarly find the largest mz peak in the doublet cluster
  cluster.max.mz <- peak.mz2
  repeat {
    higher.mz.exist <- abs(scan.mz - cluster.max.mz - isotope.mz.delta)/cluster.max.mz < mz.ppm.cut
    if ( sum( higher.mz.exist ) < 1 ) {
      # cannot find lower peak any more, finish searching
      break
    } else {
      cluster.max.mz <- scan.mz[higher.mz.exist][1]
    }
  } # repeat

  out.filename <- paste(output.path, prefix, output.file.base, "_", id,"_", peaks.n, output.file.suffix, sep="")
  # make plot for each isotopic cluster
  png(out.filename)
  par(mfrow=c(2,1))
  par(las=0)
  # mz spectrum top
  #plotScan(xfile, best.scan.number, massrange=mass.range)
  ylimit <- c( 0, max(scan.int[(scan.mz>=monoiso.mz&scan.mz<=cluster.max.mz)]) )
  plot( scan.data[,1], scan.data[,2], type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit)
  title( paste("Scan # ", best.scan.number, " @ ", round(xfile@scantime[best.scan.number], 1),
               " sec (raw # ", (best.scan.number-1)*(num.ms2.per.ms1+1)+1, " @ ",
               round(xfile@scantime[best.scan.number]/60,1)," min)", sep = "")
        )
  title( paste("Charge:",charge, "; Monoisotopic mz:", round(monoiso.mz,3),"(M)" ), line=0.5)
  #ymark <- min( scan.data[,2] ) + 10
  ymark <- 0
  #label all identified peaks under this isotopic envelope
  points( peaks.mz1, rep(ymark, length(peaks.mz1) ), pch=23, col="red",bg="red")
  points( peaks.mz2, rep(ymark, length(peaks.mz2) ), pch=23, col="blue",bg="blue")
  points( peak.mz1, peak.maxo1, pch="C", col="red",bg="red")
  points( peak.mz2, peak.maxo2, pch="C", col="blue",bg="blue")
  points( monoiso.mz, monoiso.int, pch="M", col="red",bg="red")
  #chromatogram bottom
  raw.ECI.light <- rawEIC(xfile, c(peak.mz1*(1-mz.ppm.cut), peak.mz1*(1+mz.ppm.cut)) )
  raw.ECI.heavy <- rawEIC(xfile, c(peak.mz2*(1-mz.ppm.cut), peak.mz2*(1+mz.ppm.cut)) )
  xlimit <-c( max(1, best.scan.number-25), min(best.scan.number+25, length(raw.ECI.light[[1]]) ) )
  ylimit <- range(c(raw.ECI.light[[2]][xlimit[1]:xlimit[2]], raw.ECI.heavy[[2]][xlimit[1]:xlimit[2]]))
  plot(raw.ECI.light[[1]], raw.ECI.light[[2]], type="l", col="red",xlab="scan #", ylab="intensity",
       main=paste("Chromatogram of", as.character(format(peak.mz1, digits=7)),
         "and", as.character(format(peak.mz2,digits=7)), "m/z (C)"), xlim=xlimit, ylim=ylimit)
  lines(raw.ECI.heavy[[1]], raw.ECI.heavy[[2]], col='blue', xlim=xlimit, ylim=ylimit)
  dev.off()
}

