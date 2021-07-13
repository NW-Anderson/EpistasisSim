library(data.table)
filenames <- list.files(path = "./output")

for(i in 1:length(filenames)){
  system(paste('tail -n +14 ',
               './output/',
               filenames[i],
               ' > ./trimmed/',
               substr(filenames[i],1, nchar(filenames[i]) - 3),
               'csv', sep = ''))
}

filenames <- list.files(path = './trimmed')
for(i in 1:length(filenames)){
  rawout <- as.matrix(read.csv(file = paste('./trimmed/', filenames[i], sep = '')))[,-1159]
  sigsnps6 <- vector("list",10)
  names(sigsnps6) <- paste('pop', 1:10)
  for(pop in 1:10){
    popdat <- rawout[which(rawout[,2] == pop),]
    sigsnps6[[pop]] <- which((popdat[2,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= 0.01)
  }
  jaccmat <- array(dim = c(10,10))
  for(j in 1:9){
    for(k in (j+1):10){
      jaccmat[j,k] <- length(intersect(sigsnps6[[j]], sigsnps6[[k]])) / length(union(sigsnps6[[j]], sigsnps6[[k]]))
    }
  }
  destination <- paste('./Jaccards/', substr(filenames[i],1, nchar(filenames[i]) - 3), 'RData', sep = '')
  save(jaccmat, file = destination)
}