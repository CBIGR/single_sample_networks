args <- commandArgs() #extract the argument given on the command line
print(args)
expr_table_dir <- args[6] #complete directory of input expression table
output_dir <- args[7] #directory of output table
output_name <- args[8] #name of output table
packages_dir <- args[9] #directory of R packages
lioness_dir <- args[10] #complete directory of the lioness function
netw_inf_dir <- args[11] #complete directory of the function to use for network inference


#library("GENIE3", lib.loc=packages_dir)
library("stringr", lib.loc=packages_dir)
library("data.table", lib.loc=packages_dir)
#library("tictoc", lib.loc=packages_dir)
#tic()
netw_inf_name <- str_match(netw_inf_dir, "^.*/([A-Za-z0-9_]*).R$")[,2] # !!FUNCTION MUST HAVE THE SAME NAME AS THE FILE IN WHICH THE FUNCTION IS STORED
print(netw_inf_name)

source(lioness_dir) #needed function script 1
source(netw_inf_dir) #needed function script 2
print('a')
#toc()
y_brain <- fread(file = expr_table_dir, header=TRUE, data.table=FALSE) #expression input file
print('b')
#toc()
rownames(y_brain) <- y_brain[,1]
print('c')
#toc()
y_brain[,1] <- NULL
print('d')
#toc()
y_brain <- data.matrix(y_brain)
print('e')
#toc()
y_lioness <- lioness(y_brain, f=get(netw_inf_name)) #use the lioness function implemented in the lioness.R script,
#second argument is the netFun function by default, implemented in the netFun.R script


write.csv(y_lioness, paste(output_dir, '/', output_name, netw_inf_name, '.csv', sep=""), col.names=TRUE, row.names=FALSE)
