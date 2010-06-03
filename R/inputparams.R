read.input.params <- function( file.name ) {
  params <- list()
  raw.input <- scan(file.name, what=character(), quiet=TRUE, comment.char="!")
  equal.pos <- which( raw.input == "=" )
  n.params <- length(equal.pos)
  ## param 1:(N-1)
  for ( i in 1:(n.params-1) ) {
    param.name <- raw.input[equal.pos[i]-1]
    param.data <- raw.input[(equal.pos[i]+1):(equal.pos[i+1]-2)]
    params[[param.name]] <- param.data
  }
  ## last param
  param.name <- raw.input[equal.pos[n.params]-1]
  param.data <- raw.input[(equal.pos[n.params]+1):length(raw.input)]
  params[[param.name]] <- param.data
  return(params)
}

read.chem.table <- function( table.name ) {
  orig.table <- read.table(table.name, header=T, sep="\t", comment.char="!")
  named.table <- orig.table[,2:ncol(orig.table)]
  rownames(named.table) <- orig.table[,1]
  return(named.table)
}

init.atom.mass <- function() {
  atom.mass.vec <- c(12.000000, #C
                     1.007825,  #H
                     15.994915, #O
                     14.003074, #N
                     31.972072, #S
                     30.973763, #P
                     15.000109, #N15
                     2.014102,  #H2
                     13.003355, #C13
                     1.0072765  #Hplus
                     )
  names(atom.mass.vec) <- c("C","H","O","N","S","P","N15","H2","C13","Hplus")
  return(atom.mass.vec)
}

init.aa.mass <- function(atom.mass.vec, chem.table ) {
  aa.names <- rownames(chem.table)
  chem.names <- colnames(chem.table)
  aa.mass.vec <- rep(0, length(aa.names))
  names(aa.mass.vec) <- aa.names
  for (aa in aa.names ) {
    mass <- 0.0
    for (chem in chem.names) {
      mass <- mass + atom.mass.vec[chem]*chem.table[aa,chem]
    }
    aa.mass.vec[aa] <- mass
  }
  return(aa.mass.vec)
}

calc.peptide.mass <- function(sequence, aa.mass.vec) {
  ## get rid of flanking residues deliminated by "."
  peptide.vec <- unlist( strsplit(sequence,".",fixed=T) )
  peptide.vec <- unlist( strsplit(peptide.vec[2],"",fixed=T) )
  mass <- 0
  for ( aa in peptide.vec ) {
    mass <- mass + aa.mass.vec[aa]
  }
  mass <- mass + aa.mass.vec["NTERM"] + aa.mass.vec["CTERM"]
  return(mass)
}
