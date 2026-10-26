#customGSEA, input genes as a string, e.g. c("INS", "PDX1", etc.)
customGSEA<-function(inputgenes,outputfile){
  genes<-inputgenes
  nlength<-length(genes)
  #TERM2GENE, prepare custom geneset using user input genes, geneset requires a 2 column data frame with pathway/term on the left and genes on the right 
  gene_sets<-data.frame(Term=rep("Custom Gene List", nlength), Gene=genes)
  #Load our DE lists for GSEA for each hormone into a single list of dataframes
  directory_path <- "/athena/chenlab/scratch/jjv4001/betacelldeseq"
  csv_files <- list.files(path = directory_path, pattern = "\\.csv$", full.names = TRUE)
  list_of_data_frames <- lapply(csv_files, read.csv)
  names(list_of_data_frames) <- tools::file_path_sans_ext(basename(csv_files))
  #Convert the dataframes into named vectors of decreasing log2foldchange, with genes as names
  convert_to_named_vector <- function(df, value_var, name_var) {
    vec <- setNames(df[[value_var]], df[[name_var]])
    sorted_vec <- sort(vec, decreasing = TRUE)
    return(sorted_vec)
  }
  named_vectors_list <- lapply(list_of_data_frames, convert_to_named_vector, value_var = "avg_expr", name_var = "Gene")
  order_named_vector <- function(named_vector) {
    named_vector[order(-named_vector)]
  }
  ordered_named_vectors_list <- lapply(named_vectors_list, order_named_vector)
  #Run GSEA for each hormone against custom gene list, custom gene list is supplied in TERM2GENE and genelist supplied consists of ordered named numeric vectors and create a list of the GSEA results
  run_gsea <- function(gene_list, gene_sets) {
    gsea_result <- GSEA(geneList = gene_list, TERM2GENE = gene_sets, minGSSize=5, maxGSSize=500, pvalueCutoff = 1, pAdjustMethod="BH", verbose=TRUE, nPerm=1000)
    return(gsea_result)
  }
  gsea_results_list <- lapply(ordered_named_vectors_list, run_gsea, gene_sets = gene_sets)
  for (i in seq_along(gsea_results_list)) {
    cat("GSEA Results for Vector", names(ordered_named_vectors_list)[i], ":\n")
    print(gsea_results_list[[i]])
  }
  
  filtered_gsea_results_list <- Filter(function(x) {
    dims <- dim(x)
    !(dims[1] == 0)
  }, gsea_results_list)
  #convert GSEA list to a dataframe containing GSEA results from all hormones 
  gsea_to_df <- function(gsea_result, identifier) {
    df <- as.data.frame(gsea_result)
    df$hormone <- identifier  # Add an identifier column to track the source of each result
    return(df)
  }
  gsea_dfs <- lapply(seq_along(filtered_gsea_results_list), function(i) {
    gsea_to_df(filtered_gsea_results_list[[i]], names(filtered_gsea_results_list)[i])
  })
  combined_gsea_df <- do.call(rbind, gsea_dfs)
  write.csv(combined_gsea_df,outputfile)}
  
  #Create sign column for faceting in ggplot indicating activated or suppressed if NES is above 0 or below 0, create Enrichment Score column with absolute NES and create Pathways column containing ID
  combined_gsea_df$Pathway<-combined_gsea_df$ID
  combined_gsea_df$sign<-combined_gsea_df$NES
  combined_gsea_df$`Normalized Enrichment Score`<-combined_gsea_df$NES
  combined_gsea_df$sign<-ifelse(combined_gsea_df$sign > 0, "Activated", "Suppressed")
  combined_gsea_df$EnrichmentScore<-abs(combined_gsea_df$NES)
  combined_gsea_df$`-log10(p-value)`<--log10(combined_gsea_df$pvalue)
  
  #Plot enrichment score results
  p<-ggplot(combined_gsea_df, aes(x = `Normalized Enrichment Score`, y = p.adjust, color= hormone, size=`-log10(p-value)`))+geom_point(alpha=1.0)+geom_hline(yintercept = 0.05, color = "black", linetype = "dashed", size = 0.5)+ geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p+geom_text_repel(aes(label = hormone, segment.color = NA), hjust = -0.3, vjust = 0.5, size = 4)+ guides(color = "none")

