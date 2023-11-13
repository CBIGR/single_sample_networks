args <- commandArgs() #extract the argument given on the command line
sample <- as.numeric(args[6]) #position of the wished argument
y_brain_pr <- read.csv(file = '/kyukon/scratch/gent/434/vsc43409/omics_data/JDS_brain_expression_2_75_var_only_names.csv', header=TRUE, row.names = 1)
source('/kyukon/scratch/gent/vo/000/gvo00027/vsc43409/JDS_paper_SSN/SSN_part1.R')
#print(sample)
sample_name <- colnames(y_brain_pr)
SSN_uncorrected<- SSN_edges_part1(y_brain_pr, sample)

dir <- paste('/kyukon/scratch/gent/vo/000/gvo00027/vsc43409/JDS_paper_SSN/JDS_seperate_samples/JDS_SSN_brain_2_75_uncorrected_pvalue_edges', sample_name[sample], '.txt',sep="") #file name must change per sample
write.table(SSN_uncorrected, dir, row.names = FALSE)