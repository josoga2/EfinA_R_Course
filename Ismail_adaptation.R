library(rBLAST)
library(filesstrings)

#set working directory (root directory of your files and folders)
ismailWorkingDirectory <- setwd("")

#import filenames (a txt file containing all the genome filenames as a txt file )
myGenomes <- read.delim("img_data_all_genomes/openFolders.txt", header = F)


#Metadata (if you have detailed metadata about the species)
myMetadata <- read.delim("jgi_exportdata.txt", header = T)

#move files to individual folders to create blast dbs
#do this once
for (fileN in list.files("openGenomes/")) {
  #print(fileN)
  myFiles <- strsplit(fileN, ".fna")[[1]][1]
  dir.create(paste0("FullGenomes/", myFiles))
  file.move(paste0("openGenomes/", fileN), paste0("FullGenomes/", myFiles))
}


##
#input queries (your genes)
sul1Dna <- readDNAStringSet("Queries/sul1.fasta")



#accept result for each query's percent.ident
sul1Result <- c()


#perform BLAST
for (iGen in myGenomes$V1) {
  print(which(iGen == myGenomes$V1))
  
  strainId <- c(strainId, iGen)
  
  #implement blast
  makeblastdb(paste0(suseWorkingDirectory, "/FullGenomes/", iGen, paste0("/",iGen, '.fna')))
  currBl <- blast(db = paste0(suseWorkingDirectory, "/FullGenomes/", iGen, paste0("/",iGen, '.fna')))
  
  #sul1Cl (accepts the results iteratively and incrementally)
  sul1Cl <- predict(currBl, sul1Dna)
  sul1Result <- c(sul1Result, list(sul1Cl))
  
  
}

#process the results and extract percent identity (you can extract other details too) 
sul1ResultProc <- c()

for (i in 1:length(sul1Result)) {
  if(class(sul1Result[[i]]$Perc.Ident) == 'logical'){
    sul1ResultProc <- c(sul1ResultProc, 0)
  }else{
    sul1ResultProc <- c(sul1ResultProc, sul1Result[[i]]$Perc.Ident[1])
  }
  
  
}

length(sul1ResultProc) #should be equal to the total number of input genomes


### extra processing to export your data
finalProcessedData <- data.frame("strain_id" = strainId,
                                 "sul1BLAST" = sul1ResultProc)

View(finalProcessedData)


#export data
write.table(finalProcessedData, file = "blastJGI_Output_with_metadata.csv", quote = F, row.names = F, sep = ',')

