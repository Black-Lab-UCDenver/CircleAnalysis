library(readr)
library(dplyr)
library(RBedtools)


chromSizes<- read_tsv("/path/to/hg38ChromSizes.txt", col_names = FALSE)
chromSizes$number<- gsub("chr", "", chromSizes$X1)
chromSizes<- chromSizes[order(chromSizes$number),]

totalChroms<- sum(chromSizes$X2)
for(r in 1: nrow(chromSizes)){
  if(r==1){
    start<- 1
    stop<- paste(chromSizes[r, "X2"])
  }else{
    start<- sum(chromSizes[1:(r-1), "X2"])+1
    stop<- sum(chromSizes[1:(r), "X2"])
    
  }
  
  out<- data.frame(chromSizes[r, ], start, stop)
  
  if(exists("chromFinal")==TRUE){
    chromFinal$stop<- as.numeric(chromFinal$stop)
    chromFinal<- dplyr::bind_rows(chromFinal, out)
  }else{
    chromFinal<- out
  }
  
  rm(out, start, stop, r)
}

pDir <- list.files(path = "/path/to/allCirclesByTCGAPatient", 
                   full.names = TRUE, 
                   recursive = FALSE)


for(i in pDir){
  
  setwd("/path/to/allCirclesByTCGAPatient")
  setwd(i)
  circlesIn<- read_tsv("circlesNoCentromere.bed", col_names = FALSE)
  circlesIn$length<-circlesIn$X3-circlesIn$X2
  
  name<- gsub("/path/to/allCirclesByTCGAPatient/", "", i)
  
  
  
  for(z in 1:nrow(circlesIn)){
    
    print(circlesIn[z, ])
    lengthCircles<-as.numeric(paste(circlesIn[z, 'length']))
    randomPick<-sample(1:totalChroms, 1)
    position<- dplyr::filter(chromFinal, start<randomPick & stop>randomPick)
    
    chrom<- position[1, "X1"]
    startP<- as.numeric(paste0(position$stop))-randomPick
    stopP<-startP+lengthCircles
    
    circleOut<- data.frame(chrom, startP, stopP)
    
    overLap<- data.frame(chrom, startP, stopP)%>%from_data_frame%>%
      RBedtools(tool='intersect', options= "-wa -wb", 
                a = ., 
                b = '/path/to/hg38Centromere.bed') %>%to_data_frame()
    
    
    maxLength<- position[1, 'X2']-circleOut[1, 'stopP']
    
    
    
    
    while(maxLength<0 | nrow(overLap)>0){
      rm(randomPick, position, chrom, startP, stopP, circleOut, overLap, maxLength)
      randomPick<-sample(1:totalChroms, 1)
      position<- dplyr::filter(chromFinal, start<randomPick & stop>randomPick)
      
      chrom<- position[1, "X1"]
      startP<- as.numeric(paste0(position$stop))-randomPick
      stopP<-startP+lengthCircles
      
      circleOut<- data.frame(chrom, startP, stopP)
      
      overLap<- data.frame(chrom, startP, stopP)%>%from_data_frame%>%
        RBedtools(tool='intersect', options= "-wa -wb", 
                  a = ., 
                  b = '/path/to/hg38Centromere.bed') %>%to_data_frame()
      
      
      maxLength<- position[1, 'X2']-circleOut[1, 'stopP']
      
    }
    
    chromLength<- paste0(position[1, "X2"])
    
    out<- data.frame(circleOut, lengthCircles, chromLength, circlesIn[z, ])
    
    if(exists('outFinal')==TRUE){
      outFinal<-rbind(outFinal, out)
    }else{
      outFinal<-out
    }
    
    rm(randomPick, position, chrom, startP, stopP, circleOut, overLap, 
       maxLength, chromLength, lengthCircles, out)
  }
  
  write_tsv(outFinal, paste0(i, "/", name, ".newNoCentRandCircles.bed"), col_names = FALSE)
  rm(i, circlesIn, name, position, out, circleOut, overLap, chrom, lengthCircles, maxLength, 
     randomPick, startP, stopP, z, outFinal)
}







