## file name from input args
args <- commandArgs(trailingOnly=T)

nset <- length(args)/3
ratio.files <- args[seq(1,by=3,length=nset)]
input.cols  <- args[seq(2,by=3,length=nset)]
output.cols <- args[seq(3,by=3,length=nset)]

table <- as.list(ratio.files)

all.uniq <- NULL
sp=" "
for (i in 1:nset ) {
  tmp.table <- read.table(ratio.files[i], header=T,sep="\t",quote="")
  for (j in 1:nrow(tmp.table)) {
    if ( ! is.na(tmp.table[j,"index"]) ) {
      tmp.table[j,"ipi"] = tmp.table[j+1,"ipi"]
      tmp.table[j,"description"] = tmp.table[j+1,"description"]
      tmp.table[j,"symbol"] = tmp.table[j+1,"symbol"]
    }
  }
  table[[i]] <- tmp.table[! is.na(tmp.table[,"index"]),]

  table[[i]][,"sequence"] <- as.character( table[[i]][,"sequence"] )
  table[[i]]$uniq <-i
  for (ii in 1:length(table[[i]]$ipi) ) {
    ipi <- as.character(table[[i]][ii,"ipi"])
    description <- as.character(table[[i]][ii,"description"])
    symbol <- as.character(table[[i]][ii,"symbol"])
    sequence <- as.character(table[[i]][ii,"sequence"])
    table[[i]]$uniq[ii]<- paste(ipi,description,symbol,sequence,sep=":")
  }
  all.uniq<-c(all.uniq,table[[i]]$uniq)
}

count <- 0
link.list <- as.list( levels(as.factor(all.uniq) ) )
nuniq <- length(link.list)
out.num.matrix <- matrix(NA, nrow=nuniq,ncol=1*nset)
colnames(out.num.matrix) <- output.cols
char.names <- c("index","ipi", "description", "symbol", "sequence")
out.char.matrix <- matrix(" ",nrow=nuniq,ncol=length(char.names))
colnames(out.char.matrix) <- char.names
for (uniq in levels(as.factor(all.uniq) ) ) {
  count <- count + 1
  tmp.split <- strsplit(uniq,":")[[1]]
  out.char.matrix[count,"index"] <- as.character(count)
  out.char.matrix[count,"ipi"] <- tmp.split[1]
  out.char.matrix[count,"description"] <- tmp.split[2]
  out.char.matrix[count,"symbol"] <- tmp.split[3]
  out.char.matrix[count,"sequence"] <- tmp.split[4]

  for ( i in 1:nset ) {
    match <- table[[i]][,"uniq"] == uniq
    if ( sum(match) == 1 ) {
      ratio <- table[[i]][match,input.cols[i]]
      if (ratio > 0.0) {
        out.num.matrix[count,i] <- ratio
      } else {
        out.num.matrix[count,i] <- NA
      }
    }
  }
}

## order by ratio from last set to first set
z.order <- do.call("order",c(data.frame(out.num.matrix[,seq(1,nset)]), na.last=T))

html.table <- cbind(out.char.matrix,out.num.matrix)[z.order,]

html.table[,"index"] <- seq(1,nrow(html.table))

write.table(html.table, file=paste("compare_averaged_ratios",paste(output.cols,sep="",collapse="_"),
                          "txt",sep="."),
            quote=F, sep="\t", row.names=F,na="0.00")

library(limma)
png(paste("compare_averaged_ratios_vennDiagram",paste(output.cols,sep="",collapse="_"),"png",sep="."))
venn.out.matrix <- ! is.na(out.num.matrix)
vc <- vennCounts(venn.out.matrix)
vennDiagram(vc,main="Number of peptides with valid ratios",counts.col="red")
dev.off()

if (ncol(venn.out.matrix)==2) {
  png(paste("compare_averaged_ratios_scatterplot",paste(output.cols,sep="",collapse="_"),"png",sep="."))
  out.num.matrix[is.na(out.num.matrix)] <- 0.0
  qtt <- quantile(as.numeric(out.num.matrix),probs=seq(0,1,0.01))
  limit <- c(qtt[2],qtt[length(qtt)-1])
  plot(out.num.matrix[,1], out.num.matrix[,2],main="averaged ratios comparison",xlab=output.cols[1],ylab=output.cols[2],
       xlim=limit,ylim=limit)
  abline(0,1)

  dev.off()
}
