customGSEA <- function(inputgenes, outputfile){
  
  genes <- inputgenes
  nlength <- length(genes)
  
  # custom gene set
  gene_sets <- data.frame(
    Term = rep("Custom Gene List", nlength),
    Gene = genes
  )
  
  # load DE tables
  directory_path <- "/athena/chenlab/scratch/jjv4001/betacelldeseq"
  csv_files <- list.files(path = directory_path, pattern = "\\.csv$", full.names = TRUE)
  
  list_of_data_frames <- lapply(csv_files, read.csv)
  names(list_of_data_frames) <- tools::file_path_sans_ext(basename(csv_files))
  
  # convert dataframe → ranked named vector
  convert_to_named_vector <- function(df, value_var, name_var) {
    
    df <- df[!is.na(df[[value_var]]) & !is.na(df[[name_var]]), ]
    
    vec <- setNames(df[[value_var]], df[[name_var]])
    
    # remove duplicate gene names (required for fgsea)
    vec <- vec[!duplicated(names(vec))]
    
    sort(vec, decreasing = TRUE)
  }
  
  named_vectors_list <- lapply(
    list_of_data_frames,
    convert_to_named_vector,
    value_var = "avg_expr",
    name_var = "Gene"
  )
  
  # run GSEA
  run_gsea <- function(gene_list, gene_sets){
    
    gsea_result <- GSEA(
      geneList = gene_list,
      TERM2GENE = gene_sets,
      minGSSize = 5,
      maxGSSize = 500,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE   # removes console messages
    )
    
    return(gsea_result)
  }
  
  gsea_results_list <- lapply(named_vectors_list, run_gsea, gene_sets = gene_sets)
  
  # remove empty results
  filtered_gsea_results_list <- Filter(function(x){
    if(is.null(x)) return(FALSE)
    nrow(as.data.frame(x)) > 0
  }, gsea_results_list)
  
  # convert to dataframe
  gsea_to_df <- function(gsea_result, identifier){
    
    df <- as.data.frame(gsea_result)
    df$hormone <- identifier
    df
  }
  
  gsea_dfs <- lapply(seq_along(filtered_gsea_results_list), function(i){
    gsea_to_df(filtered_gsea_results_list[[i]], names(filtered_gsea_results_list)[i])
  })
  
  combined_gsea_df <- do.call(rbind, gsea_dfs)
  
  write.csv(combined_gsea_df, outputfile)
}