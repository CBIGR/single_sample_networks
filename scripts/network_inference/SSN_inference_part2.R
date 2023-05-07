args <- commandArgs() #extract the argument given on the command line
print(args)
input_dir <- args[6] # directory folder of input files
input_pattern <- args[7] #pattern shared by all input files 
output_dir <- args[8] #directory folder and shared pattern for output files 
packages_dir <- args[9] #directory of R packages 


library("stringr", lib.loc = packages_dir)
library("data.table", lib.loc = packages_dir)
file_list <- list.files(path = paste( input_dir, '/', sep="") ,pattern = input_pattern) #all saved files with
print(file_list)

#making an empty data frame
SSN_table <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("gene_1", "gene_2", "p_value", "PCC_value")
colnames(SSN_table) <- x

#empty order vector 
order<- c()
p_values <- c()


##loop over all files stored in the previous step (for each sample 1):
for (i in 1:length(file_list)){
  print('for loop 1')
  dir <- paste( input_dir, '/', file_list[i] ,sep="") #compose directory
  ########## CHANGED (compared to SSN_inf_2_phanpy_alt.R and SSN_inf_2_phanpy_alt_samples_no_NA) ##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  sample <- str_match(file_list[i], paste( input_pattern,'(X[0-9]+).txt', sep=""))[2] # distract the sample number from the file name #####################
  ##########
  order <- c(order, sample)   #it is important to know in wich order the files (and thus samples) are looped through
  ###########################################
  table <- read.table(file = dir, header=TRUE)
  ##########################################
  print(dim(table))
  p_values <- c(p_values, table[,3]) #all p values stored in one big vector  
  #file.remove(file_list[i])
}
print(order)

print('p-adjustment') #for debugging/control on server
##the actual p-adjustment of all samples together 
p_values_adj <- p.adjust(p_values, method="BH") #p-value adjustment

##making and storing the files with corrected p values for all samples 
print(length(p_values_adj)) #for debugging/control on server

edges_per_file <- length(p_values_adj)/length(file_list) #how many edges per file/sample are present
print(edges_per_file) #for debugging/control on server

start <- 1 #a 'counter'
for (i in 1:length(order)){
  print('for loop 2')
  p_values_sample <- p_values_adj[start:(start+(edges_per_file-1))] #the subpart of the total table belonging to a sample 
  print(length(p_values_sample))#for debugging/control on server
  table_dir <- paste(input_dir, '/',input_pattern, order[i] , '.txt',sep="") #####################################"
  ###########################
  table <- read.table(file = table_dir, header=TRUE)
  ################################
  table[,3] <- p_values_sample
  start <- start+edges_per_file #counter keeping track
  print(start)
  dir_2 <- paste( output_dir, order[i], '.txt',sep="") #the wanted file 
  #name and directory 
  write.table(table, dir_2, row.names = FALSE)
}