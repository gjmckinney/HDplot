library(vcfR)
library(ggplot2)
library(dplyr)
library(stringr)

#HDplot function
HDplot<-function(vcfData){
  #set up results table
  HDplotTable<-as.data.frame(matrix(NA,nrow=dim(vcfData@gt)[1],ncol=13))
  colnames(HDplotTable)<-c("CHROM","POS","ID","depth_a","depth_b","ratio","num_hets","num_samples","num_called","H_all","H","std","D")
  
  #get genotypes from vcf file
  genos<-extract.gt(vcfData, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  
  #replace NA genotypes with ./.
  genos[is.na(genos)]<-"./."
  
  #get allele reads from vcf file
  reads<-extract.gt(vcfData, element = "AD", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
    
  #replace reads for samples with missing data with 0
  #reads<-gsub("\\.,\\.","0,0",reads)
  reads[grepl("\\.",reads)]<-"0,0"
  reads[is.na(reads)]<-"0,0"
  
  alleleReads<-apply(reads,2,function(x) str_split_fixed(x,",",2))
  alleleReads_1<-alleleReads[1:dim(reads)[1],]
  alleleReads_2<-alleleReads[dim(reads)[1]+1:dim(reads)[1],]
  #convert to numeric format
  alleleReads_1<-apply(alleleReads_1,2, function(x) as.numeric(x))
  alleleReads_2<-apply(alleleReads_2,2, function(x) as.numeric(x))
  #subset to heterozygous genotypes
  #make genotype matrix where heterozygotes are 1 and other genotypes are 0
  hetMatrix<-apply(genos,2,function(x) dplyr::recode(x,'0/0'=0,'1/1'=0,'./.'=0,'0/1'=1,'1/0'=1))
  calledGenos<-apply(genos,2,function(x) dplyr::recode(x,'0/0'=1,'1/1'=1,'0/1'=1,'1/0'=1,.default=NA_real_))
  #multiply read count matrices by heterozygote matrix to get read counts for heterozygous genotypes
  alleleReads_1_het<-alleleReads_1*hetMatrix
  alleleReads_2_het<-alleleReads_2*hetMatrix
  #rows are loci and columns are samples
  #sum reads per allele per locus for heterozygous samples
  A_reads<-apply(alleleReads_1_het,1,sum)
  B_reads<-apply(alleleReads_2_het,1,sum)
  totalReads<-A_reads+B_reads
  ratio<-A_reads/totalReads
  std<-sqrt(totalReads*0.5*0.5)
  z<- -(totalReads/2-A_reads)/std
  #get percent heterozygosity for each locus
  numHets<-apply(hetMatrix,1,sum)
  hetPerc<-numHets/dim(hetMatrix)[2]
  
  numGenos<-apply(calledGenos,1,sum,na.rm=TRUE)
  H<-numHets/numGenos
  
  #assign results to HDplotTable
  HDplotTable$CHROM<-vcfData@fix[,"CHROM"]
  HDplotTable$POS<-vcfData@fix[,"POS"]
  HDplotTable$ID<-vcfData@fix[,"ID"]
  HDplotTable$depth_a<-A_reads
  HDplotTable$depth_b<-B_reads
  HDplotTable$ratio<-ratio
  HDplotTable$num_hets<-numHets
  HDplotTable$num_samples<-dim(hetMatrix)[2]
  HDplotTable$num_called<-numGenos
  HDplotTable$H_all<-hetPerc
  HDplotTable$H<-H
  HDplotTable$std<-std
  HDplotTable$D<-z

  return(HDplotTable)
}


