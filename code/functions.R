

##-functions

##- 1. Preprocess dataframe
# This function preprocesses the input DataFrame by removing the suffix from transcript names and converting tissue names to lowercase.
preprocess_dataframe <- function(df) {
  df$transcript <- gsub("\\..*$", "", df$transcript)
  df$tissue <- tolower(df$tissue)
  return(df)
}

##-2. outersect
# This function returns elements that are in one of the sets but not both.
outersect <- function(x, y) {
  unique(c(setdiff(x, y), setdiff(y, x)))
}

##- 3.classify_promoter_events
# This function classifies promoter events as "single" or "multiple" based on the TSS differences within a gene.
classify_promoter_events <- function(df_test, threshold=100) {
  
  df = df_test %>%
    group_by(ensembl_gene_id) %>%
    arrange(TSS) %>%
    mutate(TSS_difference = abs(dplyr::first(TSS) - TSS))
  
  promoter_summary <- df %>%
    group_by(ensembl_gene_id) %>%
    summarise(promoter_event = ifelse(
      (sum(TSS_difference == 0) == 1) & all(TSS_difference[TSS_difference != 0] < threshold),
      "single",
      "multiple"
    ))
  
  promoter_summary_df = merge(promoter_summary, df, by = "ensembl_gene_id")
  promoter_summary_df <- promoter_summary_df %>% as.data.frame %>% dplyr::select(-promoter_event, everything(), promoter_event)
  
  return(promoter_summary_df)
}


##- 4. calculate_AFE_ALE 
# This function calculates AFE (Alternative First Exon) and ALE (Alternative Last Exon) events.
calculate_AFE_ALE <- function(df) {
  
  df_filtered <- df %>%
    # filter(status == "exon") %>%
    mutate(start = as.integer(TranscriptStart), end = as.integer(TranscriptEnd)) %>%  
    group_by(ensembl_gene_id, strand) %>%
    mutate(
      AFE_event = if_else(any(TranscriptStart != dplyr::first(TranscriptStart)), "yes", "no"),  
      ALE_event = if_else(any(TranscriptEnd != dplyr::first(TranscriptEnd)), "yes", "no")  
    ) 
  return(df_filtered)
}


##- 5. count_AFE_ALE_genes
# This function counts the number of genes and transcripts with AFE or ALE events.
count_AFE_ALE_genes <- function(df, col1="ALE_event") {
  
  event_genes <- df %>% dplyr::filter(!!sym(col1) == "yes")
  
  gene_n <- event_genes$ensembl_gene_id %>% unique() %>% length()
  tc_n <- event_genes$transcript %>% unique() %>% length()
  
  print(paste("Number of transcripts with an", col1, "event:"))
  print(tc_n)
  
  print(paste("Number of genes with an", col1, "event:"))
  print(gene_n)
  
}


##- 6. exon skip events
# This function counts exon skipping events.
count_ES_events <- function(df, 
                            strand_col = "strand", 
                            start_col = "start", 
                            end_col = "end", 
                            status_col = "status", 
                            transcript_col = "transcript") {
  
  # Identify unique exons across all transcripts for each gene
  unique_exons <- df %>% 
    filter(!!sym(status_col) == 'exon') %>%
    group_by(ensembl_gene_id, !!sym(start_col), !!sym(end_col)) %>%
    summarise(count = n(), transcripts = list(unique(!!sym(transcript_col))), .groups = 'drop') %>%
    filter(count == 1)
  
  # Identify introns
  introns <- df %>% 
    filter(!!sym(status_col) == 'intron') %>%
    dplyr::select(ensembl_gene_id, !!sym(start_col), !!sym(end_col))
  
  # Remove overlapping exons with introns
  joined_df <- full_join(unique_exons, introns, by = "ensembl_gene_id", suffix = c("_exon", "_intron"))
  overlap_condition <- joined_df$start_exon <= joined_df$end_intron & joined_df$end_exon >= joined_df$start_intron
  non_overlapping_exons <- joined_df %>%
    filter(!overlap_condition) %>% 
    dplyr::select(ensembl_gene_id, start_exon, end_exon, count, transcripts) %>%
    distinct()
  
  # Count ES events for each gene
  ES_counts <- non_overlapping_exons %>%
    group_by(ensembl_gene_id) %>%
    summarise(ES_count = n(), .groups = 'drop')
  
  return(ES_counts)
}



##- 7. intron retention
# This function identifies intron retention events by checking overlaps with exons.
has_IR_events <- function(df) {
  
  # Separate into exons and introns
  df_exons <- df %>% filter(status == 'exon')
  df_introns <- df %>% filter(status == 'intron')
  
  # Identify intron retention by checking overlaps
  IR_events <- df_introns %>%
    rowwise() %>%
    mutate(has_IR = any((start >= df_exons$start & start <= df_exons$end) &
                          (end >= df_exons$start & end <= df_exons$end))) %>%
    ungroup()
  
  # Summarize for each transcript
  IR_classification <- IR_events %>%
    group_by(ensembl_gene_id, transcript) %>%
    summarise(has_IR = any(has_IR)) %>% 
    ungroup()
  
  return(IR_classification)
}

# Longer version to count intron retention events

count_IR_events <- function(df) {
  
  # Separate into exons and introns
  df_exons <- df %>% filter(status == 'exon')
  df_introns <- df %>% filter(status == 'intron')
  
  # Cross join introns and exons, then filter and count
  IR_events <- df_introns %>%
    inner_join(df_exons, by = "ensembl_gene_id", suffix = c("_intron", "_exon")) %>%
    filter(start_intron >= start_exon & end_intron <= end_exon) %>%
    group_by(ensembl_gene_id) %>%
    summarise(IR_count = n_distinct(start_intron, end_intron)) 
  
  return(IR_events)
}


##- 8. alternate donor and acceptor splice
# This function counts alternative 3' and 5' splice sites.
count_alt_3_5_splice_sites <- function(df) {
  
  # Copy the input dataframe
  df_exons = df
  
  # Identify alternative 5' splice sites (alt_5ss)
  # Group the data by gene ID and start coordinate, then arrange by start position
  # If there are more than one distinct end positions for the same start position, label as "Yes"
  df_exons <- df_exons %>%
    group_by(ensembl_gene_id, start) %>%
    arrange(start) %>%
    mutate(
      alt_5ss = if_else(n_distinct(end) > 1, "Yes", "No")
    ) %>%
    ungroup() # Remove the grouping
  
  # Identify alternative 3' splice sites (alt_3ss)
  # Group the data by gene ID and end coordinate
  # If there are more than one distinct start positions for the same end position, label as "Yes"
  df_exons <- df_exons %>%
    group_by(ensembl_gene_id, end) %>%
    mutate(
      alt_3ss = if_else(n_distinct(start) > 1, "Yes", "No")
    ) %>%
    ungroup() # Remove the grouping
  
  # Return the modified dataframe
  return(df_exons)
}


##- 9. exon intron plot function
# This function plots the gene structure of exons and introns for a given gene.
plot_gene_structure <- function(df, gene_id, status_col = "status", start_col = "start", end_col = "end", transcript_col = "transcript") {
  # Filter the DataFrame for the specified gene ID
  exampledf <- df[df$ensembl_gene_id == gene_id, ]
  
  # Filter out exons and introns
  exons <- exampledf %>% dplyr::filter(!!sym(status_col) == "exon")
  introns <- exampledf %>% dplyr::filter(!!sym(status_col) == "intron")
  
  # Create the plot
  p <- ggplot() +
    geom_rect(data = exons, aes(
      xmin = !!sym(start_col),
      xmax = !!sym(end_col),
      ymin = as.numeric(factor(!!sym(transcript_col))) - 0.4,
      ymax = as.numeric(factor(!!sym(transcript_col))) + 0.4,
      fill = !!sym(status_col)
    )) +
    geom_segment(data = introns, aes(
      x = !!sym(start_col),
      xend = !!sym(end_col),
      y = as.numeric(factor(!!sym(transcript_col))),
      yend = as.numeric(factor(!!sym(transcript_col)))
    )) +
    scale_y_discrete(name = "Transcripts") +
    scale_x_continuous(name = "Genomic Coordinates") +
    ggtitle(paste("Gene ID:", unique(exampledf$ensembl_gene_id))) +
    theme_minimal()
  
  return(p)
}



