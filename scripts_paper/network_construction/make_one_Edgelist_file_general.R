#' 
#' makes one file of seperate edgelist files 
#' 
#' 
#' @param dir the directory of the files 
#' @param pattern the pattern the files have in common
#' @param file_ext the extension of the files with the given input pattern (e.g txt, csv)
#' @param file_name which file name to give to the made file (without extension)
#' @param col_n which column in the files with given input pattern are the wanted edge weights
#' @param p_values if TRUE there are also p_values present in the files 
#' @param p_cutoff the wanted cutoff on those p_value, all weights belonging to a p_values bigger than this cutoff will be made equal to zero 
#' @param col_p_n the column where the p_values can be found in the files 
#' @param n_genes the amount of genes in the files
#' 
make_one_edgelist_file_general <- function(dir, pattern, file_ext, file_name, col_n, p_values=FALSE, p_cutoff=0.05, col_p_n=0, n_genes){
  if (substr(dir, nchar(dir), nchar(dir)) != '/'){
    dir <- paste(dir, '/', sep="")
  }
  print(dir)
  file_list <- list.files(path = dir ,pattern = pattern)
  print(file_list)
  
  table_total <- data.frame(matrix(ncol=length(file_list)+2, nrow=n_genes^2))
  col_names <- c("reg", "tar")
  for (i in 1:length(file_list)){
    
    sub_in <- paste(pattern, '([0-9]+)', '.', file_ext)
    sample <-  sub(sub_in, "\\1", file_list[i])
    print(file_list[i]) 
    print(paste(dir, file_list[i], sep=""))
    table_sample <- fread(paste(dir, file_list[i], sep=""), header=TRUE, data.table=FALSE) ####
    #table_sample <- data.frame(table_sample) ####
    print(dim(table_sample))
    
    if(i==1){
      table_total[,1] <- table_sample[,1]
      table_total[,2] <- table_sample[,2]
      print('A')
    }
    if(i!=1){
      if(table_total[,1] != table_sample[,1] || table_total[,2] != table_sample[,2]){
        stop('there are files that do not have the edges in same order')
        print('B')
      }
    }
    
    if(p_values==TRUE){
      table_sample[which(table_sample[,col_p_n]>p_cutoff),col_n] <- 0
    }
    print('C')
    sample_name <- str_match(colnames(table_sample)[col_n], "^.*(X[0-9]+).*$")[2] ### CHANGED
    col_names <- c(col_names, sample_name)
    print(sample_name)
    print(col_names)
    print('E')
    print(match(sample_name, col_names))
    table_total[,match(sample_name, col_names)] <- table_sample[,col_n]
    
  }
  colnames(table_total) <- col_names
  fwrite(table_total, file=paste(dir, file_name, '.csv', sep=""), row.names=FALSE)
}