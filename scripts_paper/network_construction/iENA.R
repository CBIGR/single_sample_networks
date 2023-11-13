args <- commandArgs() # the input given in the command line on the hpc server when running the jobscript 
sample <- as.numeric(args[6]) # getting the right input argument represeting the sample number 
y_brain_pr <- read.csv(file = "/kyukon/scratch/gent/434/vsc43409/omics_data/JDS_brain_expression_2_75_var_only_names.csv", header=TRUE, row.names = 1)
source('/scratch/gent/vo/000/gvo00027/vsc43409/JDS_paper_iENA/iENA.R') #the function script 

sPCC_edges<- iENA(y_brain_pr, sample) #using the implemented function 

dir <- paste('/scratch/gent/vo/000/gvo00027/vsc43409/JDS_paper_iENA/iENA_brain_2_75_sample_', sample, '.csv',sep="") #output directory 
write.csv(sPCC_edges, dir, row.names = FALSE) 