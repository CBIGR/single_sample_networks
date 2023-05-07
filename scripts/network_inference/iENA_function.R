#' 
#' Calculates the sPCC as proposed in the iENA paper 
#' 
#' @param y the dataframe with preprocessed expression data; the rows are genes and the columns samples 
#' @param sample for which sample (number according to the columns in y) must the sPCC be calculated 

iENA <- function(y, sample){
  sample_name <- colnames(y)[sample]
  print(sample_name)
  y_m <- data.matrix(y)
  mean <- apply(y_m, 1, mean) # get means for all genes (rows)
  variance <- apply(y_m, 1, var) # get variances for all genes (rows)
  y_char <- data.frame(mean, variance)
  
  # THe following implementation is the sPCC formula as proposed in the iENA paper. 
  # All the input samples themselves are taken as reference samples. 
  n <- dim(y)[1] #amount of genes 
  Xi <- matrix(y_m[,sample], n, n) #expression of gene i, same gene within a row (each row containss all the same values)
  Xj <- matrix(rep(y_m[,sample],each=n),nrow=n) #expression of gene j, same gene within a column (each column contains all the same values)
  Mi <- matrix(y_char$mean, n, n) #average of gene i, same gene within a row
  Mj <- matrix(rep(y_char$mean,each=n),nrow=n) #average of gene j, same gene within a column
  V1 <- matrix(y_char$var, n, n) #variance gene i, same gene within a row
  V2 <- matrix(rep(y_char$var,each=n),nrow=n) #variance gene j, same gene within a column
  #the actual formula
  sPCC <- (Xi-Mi)*(Xj-Mj)/sqrt(V1*V2) 
  
  rownames(sPCC) <- rownames(y)
  colnames(sPCC) <- rownames(y)
  
  #making an edgelist as output
  c1 <- rep(colnames(sPCC), length(colnames(sPCC)))
  c2 <- rep(rownames(sPCC), each = length(rownames(sPCC)))
  sPCC_edges <- data.frame(c1, c2, c(sPCC))
  colnames(sPCC_edges) <- c('reg', 'tar', sample_name)
  
  return(sPCC_edges)
}