# Written by Dr. A M Mahedi Hasan, Post-doctoral Research Associate, The Leach Lab, ICB, University of Edinburgh, UK.
# This script is under GPLv3 licensing criteria.

# This script deals with comparative ChIP-seq data analysis with a combination of library size normalisation (suggested by Simon Anders & Huber, 2010) and further normalisation of IP data by input data (where available).
# In addition, R/Bioconductor package DESeq is used here for the comparative analysis of IP data (with biological replicates) between different conditions or different strains of E. coli. 
# For the later case, at least three biological replicates of each strain or condition are needed for a statistically significant fold-change measurement of DNA enrichment at a defined region in the genome, otherwise, the same method can be used for exploratory data analysis without any statistical significance. In the latter case, a minimum of two independent biological repeats are recommended to confirm qualitative reproducibility.


require(ggplot2)
require(seqinr)
require(DESeq)



# Loading the data
data <- read.csv("1kb_binned_count_data.csv")


# A library size normalisation method, suggested by Simon Anders & Huber, 2010
median_normalization <- function(data){
  data_gm <- transform(data, gm=apply(data, 1, function(data){exp(mean(log(data), na.rm = T))}))
  data_div_by_gm <- apply(data, 2, function(data){data/data_gm$gm})
  sizing_factor <- apply(data_div_by_gm, 2, median,na.rm=T)
  norm_data <- data.frame(t(apply(data, 1, function(data){data/sizing_factor})))
  return(norm_data)
}


normalised_data <- median_normalization(data[,2:9])

plot_data <- data.frame(position=seq(from= 500,by = 1000, length.out = nrow(data)))

plot_data <- cbind(plot_data, normalised_data)

plot_data <- transform(plot_data, DL4184A=plot_data$DL4184A_IP/plot_data$DL4184A_IN)
plot_data <- transform(plot_data, DL4184B=plot_data$DL4184B_IP/plot_data$DL4184B_IN)
plot_data <- transform(plot_data, DL4201A=plot_data$DL4201A_IP/plot_data$DL4201A_IN)
plot_data <- transform(plot_data, DL4201B=plot_data$DL4201B_IP/plot_data$DL4201B_IN)

plot_data <- transform(plot_data, avg_DL4184IP=(plot_data$DL4184A_IP+plot_data$DL4184B_IP)/2)
plot_data <- transform(plot_data, avg_DL4201IP=(plot_data$DL4201A_IP+plot_data$DL4201B_IP)/2)

plot_data <- transform(plot_data, avg_DL4184=(plot_data$DL4184A+plot_data$DL4184B)/2)
plot_data <- transform(plot_data, avg_DL4201=(plot_data$DL4201A+plot_data$DL4201B)/2)

plot_data <- transform(plot_data, avg_DL4184_vs_avg_DL4201=plot_data$avg_DL4184/plot_data$avg_DL4201)

plot_data <- transform(plot_data, avg_DL4184_IP=(plot_data$DL4184A_IP+plot_data$DL4184B_IP)/2)
plot_data <- transform(plot_data, avg_DL4201_IP=(plot_data$DL4201A_IP+plot_data$DL4201B_IP)/2)

plot_data <- transform(plot_data, avg_IP_DL4184_vs_avg_IP_DL4201=plot_data$avg_DL4184_IP/plot_data$avg_DL4201_IP)


##############################
#######DESeq#######################

# A: DL4184IP vs DL4201IP

countData_A <- data.frame(data[,c(3,5,7,9)], row.names = data$Bin)
condition <- factor(c("DL4184","DL4184","DL4201","DL4201"))
DDS <- newCountDataSet(countData_A, condition)
DDS <- estimateSizeFactors(DDS)
sizeFactors(DDS)
DDS2 <- estimateDispersions(DDS,fitType="local")
res <- nbinomTest(DDS2,"DL4201","DL4184")


plot_data <- cbind(plot_data, DESeq_DL4184IP_vs_DL4201IP=res$foldChange)
plot_data <- transform(plot_data, ideal_vs_DESeqIPcomp = plot_data$avg_DL4184_vs_avg_DL4201/plot_data$DESeq_DL4184IP_vs_DL4201IP)
# B: DL4184IP vs DL4184IN

countData_B <- data.frame(data[,c(3,5,2,4)], row.names = data$Bin)
condition <- factor(c("DL4184IP","DL4184IP","DL4184IN","DL4184IN"))
DDS <- newCountDataSet(countData_B, condition)
DDS <- estimateSizeFactors(DDS)
sizeFactors(DDS)
DDS2 <- estimateDispersions(DDS,fitType="local")
res <- nbinomTest(DDS2,"DL4184IN","DL4184IP")


plot_data <- cbind(plot_data, DESeq_DL4184IP_vs_DL4184IN=res$foldChange)


# C: DL4201IP vs DL4201IN

countData_C <- data.frame(data[,c(7,9,6,8)], row.names = data$Bin)
condition <- factor(c("DL4201IP","DL4201IP","DL4201IN","DL4201IN"))
DDS <- newCountDataSet(countData_C, condition)
DDS <- estimateSizeFactors(DDS)
sizeFactors(DDS)
DDS2 <- estimateDispersions(DDS,fitType="local")
res <- nbinomTest(DDS2,"DL4201IN","DL4201IP")


plot_data <- cbind(plot_data, DESeq_DL4201IP_vs_DL4201IN=res$foldChange)





