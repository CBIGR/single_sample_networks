SSN_edges_part1 <- function(expr_data, sample){
  sample_name <- colnames(expr_data)[sample] ###
  print(sample_name) ###
  expr_data_m <- data.matrix(expr_data)
  print(dim(expr_data_m))
  n <- dim(expr_data_m)[2]
  PCC_ref <- cor(t(expr_data_m)) #pearson correlation network of all samples
  PCC_ref_zs <- cor(t(expr_data_m[,-sample])) # pearson correlation network of all samples except the sample of interest 
  deltaPCC_s <- PCC_ref - PCC_ref_zs #difference in pearson correlation 
  Zstat_s <- deltaPCC_s/(1-(PCC_ref_zs^2))*(n-2) #-2 because n = all samples -1; z_statistic as defined in the Liu et al. paper (2016)
  print(dim(Zstat_s))
  Zstat_s_vect <- c(Zstat_s) #making a vector 
  print(dim(Zstat_s_vect))
  pstat_s_o <- 2*pnorm(-abs(Zstat_s_vect)) #-abs() because when the Zscore is positive the actual pvalue is 1-foundpvalue 
  pstat_s_o[pstat_s_o>1] <- 1 #all values bigger than one should be set equal to one 
  
  #getting the right rownames for the vector that will be formed by the p.adjust function 
  rownames_1 <- rep(row.names(Zstat_s), length(row.names(Zstat_s)))
  rownames_2 <- rep(row.names(Zstat_s), each=length(row.names(Zstat_s)) )
  #rownames_p <- paste(rownames_1, rownames_2, sep="_")
  pstat_s <- data.frame(rownames_1, rownames_2,pstat_s_o, c(deltaPCC_s)) #collect all info of an edge in one data frame
  #print(pstat_s)
  colnames(pstat_s) <- c("gene_1", "gene_2", "p_value", paste("PCC_value", sample_name, sep=" ")) ###
  return(pstat_s)
}