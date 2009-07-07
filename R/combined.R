uniq.tryptic.sequence <- function( raw.sequence ) {
  ## find the smallest tryptic fragment with the labeling sites
  seq.vec <- unlist(strsplit(raw.sequence,""))
  label.pos <- which(seq.vec=="*")
  if ( length(label.pos) == 0 ) {
    ## no mod site, return center sequence
    return( strsplit(raw.sequence,".",fixed=T)[[1]][2] )
  } else {
    first.label <- min( label.pos )
    start <- first.label - 1
    while(seq.vec[start] != ".") {
      start <- start - 1
      if ( (seq.vec[start] == "R" | seq.vec[start] == "K") & seq.vec[start+1] != 'P' ) {
        break
      }
    }
    last.label <- max( label.pos )
    end <- last.label+1
    while( seq.vec[end] != "." ) {
      end <- end + 1
      if ( (seq.vec[end-1] == "R" | seq.vec[end-1] == "K") & seq.vec[end] != 'P' ) {
        break
      }
    }
    uniq.seq.vec <- seq.vec[(start+1):(end-1)]
    return( paste( uniq.seq.vec[uniq.seq.vec != "*"],sep="",collapse="") )
  }
}

## file name from input args
args <- commandArgs(trailingOnly=T)

input.file <- args[1]
dirs <- args[-1]
table <- as.list(dirs)
r2.cutoff <- 0.8

## rename columns
v1 <- c(9,10,11)
vn1 <- c("IR.1.1",  "IR.1.5", "IR.1.10" )
v2 <- c(12,13,14)
vn2 <- c("LR.1.1",  "LR.1.5", "LR.1.10" )
v3 <- c(15,16,17)
vn3 <- c("R2.1.1",  "R2.1.5", "R2.1.10" )

all.table <- NULL
for (i in 1:length(dirs) ) {
  table[[i]] <- read.table(paste(dirs[i],input.file,sep=""),header=T,sep="\t",quote="")
  names(table[[i]])[c(v1,v2,v3)] <- c(vn1,vn2,vn3)
  table[[i]][,"sequence"] <- as.character( table[[i]][,"sequence"] )
  table[[i]]$run<-i
  ##table[[i]]$ipi<-i
  table[[i]]$uniq <-i
  table[[i]]$filter <- 0
  for (ii in 1:length(table[[i]]$ipi) ) {
    ##table[[i]]$ipi[ii] <- strsplit(as.character(table[[i]][ii,1]), " - ")[[1]][1]
    sequence <- as.character(table[[i]][ii,"sequence"])
    ##sequence <- strsplit(sequence,".",fixed=T)[[1]][2]
    ##sequence <- strsplit(sequence,"*",fixed=T)[[1]]
    ##sequence <- paste(sequence,sep="",collapse="")
    ##sequence <- strsplit(sequence,"")[[1]]
    ##if ( sequence[1] == "R" | sequence[1] == "K" ) {
    ##  sequence <- sequence[-1]
    ##}
    ##sequence <- paste(sequence,sep="",collapse="")
    sequence <- uniq.tryptic.sequence( sequence )
    table[[i]]$uniq[ii]<- sequence
  }
  all.table<-rbind(all.table, table[[i]])
}

## set to NA if not passing r2.cutoff
for( i in 1:length(vn1) ) {
  invalid <- (all.table[[vn3[i]]]<r2.cutoff)
  all.table[[vn1[i]]][invalid] <- NA
  #invalid <- (all.table[[v3[i]]]<=0.0)
  #all.table[[v3[i]]][invalid] <- NA
}
## only consider entries with at least two valid ratios out of three concentrations
for( i in 1:nrow(all.table) ) {
  all.table[i,"filter"] <- sum(all.table[i,vn1]>0, na.rm=T)
}

out.table <- data.frame(Count=character(0), IPI=character(0), Description=character(0),
                        Symbol=character(0), Sequence=character(0),
                        ##cr.1.1=numeric(0), cr.1.5=numeric(0), cr1.10=numeric(0),
                        mass=character(0),
                        mr.1.1=numeric(0), mr.1.5=numeric(0), mr.1.10=numeric(0),
                        sd.1.1=numeric(0), sd.1.5=numeric(0), sd.1.10=numeric(0),
                        run=numeric(0), charge=numeric(0), segment=numeric(0),
                        link=character(0)
                        )
sp=" "
count <- 0
link.list <- as.list( levels(as.factor(all.table$uniq) ) )
for (uniq in levels(as.factor(all.table$uniq) ) ) {
  ipi <- strsplit(uniq,":")[[1]][1]
  seq <- strsplit(uniq,":")[[1]][2]
  match <- all.table[,"uniq"] == uniq  ##(all.table[,"sequence"]==seq) & (all.table$ipi==ipi)
  sub.table <- all.table[match,]
  s1 <- sub.table[,"ipi"]
  s2 <- sub.table[,"sequence"]
  s5 <- -sub.table[,"filter"]
  s3 <- sub.table[,"charge"]
  s4 <- sub.table[,"segment"]

  ii <- order(s1,s2,s5,s3,s4)
  count <- count+1
  link.list[[count]] <- which(match)[ii]

  pass <- sub.table$filter>=2
  c0 <- as.character(count)
  c1 <- sp ##as.character(unique(sub.table[,"ipi"])[1])
  c2 <- sp ##as.character(unique(sub.table[,"description"])[1])
  c3 <- sp ##as.character(unique(sub.table[,"symbol"])[1])
  c4 <- as.character(uniq)
  ##c4 <- round(mean(sub.table[,"cr.1.1"],na.rm=T),digits=2)
  ##c5 <- round(mean(sub.table[,"cr.1.5"],na.rm=T),digits=2)
  ##c6 <- round(mean(sub.table[,"cr.1.10"],na.rm=T),digits=2)
  c7 <- sp ##round(mean(sub.table[,"mass"],na.rm=T),digits=2)
  c8 <- round(median(sub.table[pass,vn1[1]],na.rm=T),digits=2)
  c9 <- round(median(sub.table[pass,vn1[2]],na.rm=T),digits=2)
  c10 <- round(median(sub.table[pass,vn1[3]],na.rm=T),digits=2)
  c11 <- round(sd(sub.table[pass,vn1[1]],na.rm=T),digits=2)
  c12 <- round(sd(sub.table[pass,vn1[2]],na.rm=T),digits=2)
  c13 <- round(sd(sub.table[pass,vn1[3]],na.rm=T),digits=2)
  c14 <- as.numeric(paste(levels(as.factor(sub.table[,"run"])),sep="",collapse=""))
  c15 <- as.numeric(paste(levels(as.factor(sub.table[,"charge"])),sep="",collapse=""))
  c16 <- as.numeric(paste(levels(as.factor(sub.table[,"segment"])),sep="",collapse=""))
  this.entry <- data.frame(Count=c0,IPI=c1, Description=c2, Symbol=c3, Sequence=c4,
                           ##cr.1.1=c4, cr.1.5=c5, cr1.10=c6,
                           mass=c7,
                           mr.1.1=c8, mr.1.5=c9, mr.1.10=c10,
                           sd.1.1=c11, sd.1.5=c12, sd.1.10=c13,
                           run=c14, charge=c15,segment=c16,link=sp
                           )
  out.table <- rbind(out.table,this.entry)
}

s1 <- out.table[,"mr.1.10"]
s2 <- out.table[,"mr.1.5"]
s3 <- out.table[,"mr.1.1"]
s4 <- out.table[,"mass"]
ii <- order(s1, s2, s3, s4)
out.table <- out.table[ii,]

html.table <- out.table[F,]
uniq.count <- 0
for ( m in 1:nrow(out.table) ) {
  index <- as.numeric(out.table[m,"Count"])
  match <- link.list[[index]]
  sub.table <- all.table[match,]

  links <- sub.table[,"link"]
  runs <- as.numeric(sub.table[,"run"])

  this.entry <- out.table[m,]
  html.table <- rbind(html.table, this.entry)
  ## count only unique peptide with unique labelling site
  sequence.match <- which(as.character(this.entry[,"Sequence"]) == as.character(html.table[,"Sequence"]))
  if ( length(sequence.match) == 1 ) {
    this.count <- uniq.count <- uniq.count+1
  } else {
    this.count <- html.table[sequence.match[1],"Count"]
  }
  html.table[nrow(html.table),"Count"] <- as.character(this.count)
  for ( l in 1:length(links) ) {
    linkfile <- strsplit(as.character(links[l]),'"')[[1]]
    new.filename <- paste(dirs[runs[l]],linkfile[2],sep="")
    new.count <- paste(this.count, l, sep=".")
    new.link <- paste('=HYPERLINK(\"',new.filename,'\",\"',new.count,'\")',sep='')
    ##new.count <- paste(m, l, sep=".")
    this.entry <- data.frame(Count=sp,
                             IPI=as.character(sub.table[l,"ipi"]),
                             Description=as.character(sub.table[l,"description"]),
                             Symbol=as.character(sub.table[l,"symbol"]),
                             Sequence=as.character(sub.table[l,"sequence"]),
                             ##cr.1.1=sub.table[l,"cr.1.1"],
                             ##cr.1.5=sub.table[l,"cr.1.5"],
                             ##r1.10=sub.table[l,"cr.1.10"],
                             mass=as.character(sub.table[l,"mass"]),
                             mr.1.1=sub.table[l,vn1[1]],
                             mr.1.5=sub.table[l,vn1[2]],
                             mr.1.10=sub.table[l,vn1[3]],
                             sd.1.1=NA, sd.1.5=NA, sd.1.10=NA,
                             run=sub.table[l,"run"],
                             charge=sub.table[l,"charge"],
                             segment=sub.table[l,"segment"],
                             link=new.link
                           )
    html.table <- rbind(html.table,this.entry)
  }
}
write.table(html.table,file="combined.txt", quote=F, sep="\t", row.names=F,na="0.00")
