args <- commandArgs() #extract the argument given on the command line
print(args)
expr_table_dir <- args[6] #complete directory of input expression table
output_dir <- args[7] #directory (+ name) of output table 
packages_dir <- args[8] #directory of R packages 
functions_dir <- args[9] #complete directory of the function file of SSPGI: found on the SSPGI GitHub page https://github.com/Marscolono/SSPGI
net_dir_or_FALSE <- args[10] # complete directory of the 'backgroun net' as described in the paper (\t seperator), or FALSE if all possible edges should be considered

########################################expression data and network data
library(stringr)
library(data.table)

Tumor_expr<- fread(expr_table_dir, data.table=F, header=T)
print(Tumor_expr[1:20,1])
rownames(Tumor_expr) <- Tumor_expr[,1] #put first column as row names
Tumor_expr <- Tumor_expr[,-1]
print('A')


if (net_dir_or_FALSE == FALSE){
  ## build background 'net' W
  genes <- rownames(Tumor_expr)
  all_genes1 <- rep(genes, each=length(genes))
  all_genes2 <- rep(genes, length(genes))
  net <- cbind(all_genes1,all_genes2)
} else {
  net <- read.table(net_dir_or_FALSE, sep= "\t", header = F)
}


########################################convert expression matrix to delta rank matrix
source(functions_dir)

Tumor_expr_mc <- cbind(Tumor_expr, apply(Tumor_expr, 1, mean)) # the matrix with a mean column added 

rank_Tumor <- rank.matrix(Tumor_expr)
rank_Tumor_mc <- rank.matrix(Tumor_expr_mc)
n = dim(Tumor_expr)[2]
x = cbind(rank_Tumor_mc, rank_Tumor)
print(dim(x))



n_ref = n
n_normal = n 
n_cancer = n
deltarank.result <- delta.rank(net, x, n_ref, n_cancer) 
#save(deltarank.result, file="./deltarank.result.Rdata")

####################################### caculate the edge-perturbation matrix
EPm_normal <- EPm (deltarank.result[[2]], deltarank.result[[3]])
EPm_cancer <- EPm (deltarank.result[[2]], deltarank.result[[4]]) 
all.equal(EPm_normal, EPm_cancer) #should be TRUE in our case ####################"

#some control statements
print(all.equal(EPm_normal, EPm_cancer))
print(dim(EPm_cancer))
dim(EPm_cancer)[1] == dim(data.frame(deltarank.result[[1]]))[1] #amount of rows = amount of edges in deltarank.result[[1]]

edges <- deltarank.result[[1]]
colnames(edges) <- c("node1", "node2")
EPm_cancer_table <- cbind(edges, EPm_cancer)
#fwrite(EPm_cancer, file=output_dir)
write.table(EPm_cancer_table, file=output_dir, row.names=TRUE, col.names=NA)


#####################################################################################################################################








