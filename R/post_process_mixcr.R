#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' post_process_mixcr
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title post_process_mixcr
#' 
#' @description 
#' Joins individual sample files into one tsv file.  Also converts ucsc names to hgnc_symbol|entrez_ids
#' and outputs the gene level counts.
#' 
#' @param input_file_paths Character vector of paths to the pipeline output data. 
#' @param output_dir Path to the output folder.
#' @param sample_data_table Data table which should contatin the sample_id_column and sample_folder_column
#' @param sample_folder_column The name of the column that has sample folder names
#' @param sample_id_column The name of the column that has sample ids
#' @param thread_num Integer number of threads to run mclapply statements
#' 
#' @return A path to the rds file.
#' 
#' @export
post_process_mixcr = function(
  input_file_paths,
  output_dir,
  sample_data_table,
  my_chains = my_chains,
  sample_folder_column = "Sample_Folder",
  sample_id_column = "Sample_ID",
  thread_num = 1
){
  library(binfotron)
  
  run_common_clone_steps = function(this_df){
    # drop incomplete sequences
    this_df =  this_df[!grepl("*", this_df$ntCDR3, fixed = T), ]
    this_df =  this_df[!grepl("_", this_df$ntCDR3, fixed = T), ]
    
    # add length
    this_df$Length = nchar(this_df$ntCDR3)
    return(this_df)
  }
  
  format_regions = function(my_region){
    asterisk_location = gregexpr(pattern ='*',my_region, fixed = T)[[1]][1]
    substring(my_region, 1, asterisk_location-1)
  }
  
  format_all_regions = function(my_regions, my_split = "_"){
    output_vector = sapply(my_regions, function(x){
      all_regions = strsplit(x, split = my_split)[[1]]
      formated_regions = unique(sapply(all_regions, format_regions))
      return(paste0(formated_regions, collapse = ","))
    })
    return(output_vector)
  }
  
  readme_path = file.path(output_dir, "readme.txt")
  if(file.exists(readme_path)){ file.remove(readme_path)}
  
  a = function(...){
    my_output = paste0(...)
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    message(my_output)
  }
  
  a("Post processing MiXCR diversity using MiXCR::post_process_mixcr v", packageVersion("MiXCR"))
  a("")
  a("Diversity metrics are calculated on each chain type seperately.")
  a("NT metrics are by nucleotide sequence")
  a("AA metrics are by amino acid sequence")
  a("VRegion metrics are by V region (top hit) diversity")
  a("JRegion metrics are by J region (top hit) diversity")
  a("Chain abundance and fraction are the same across NT, AA, VRegion, and JRegion")
  a("")
  a("Definitions")
  a("* Abundance: reads found for all cdr3 of that chain")
  a("* Richness: unique cdr3 found for that chain")
  a("* Fraction: fraction of reads for that chain")
  a("* Convergence: (number of unique ntCDR3 that share an aaCDR3 with at least one other unique ntCDR3) / (ntCDR3 Richness)")
  a("* dXX: 1 - n / N, where n is the minimum number of clonotypes accounting for at least XX% of the total reads and code N is the total number of clonotypes in the sample.")
  a("* Chao1: from vegan::estimateR 'S.chao1'")
  a("* Inv_Simpson: from vegan::diversity 'invsimpson'")
  a("* Shannon_Entropy: from vegan::diversity 'shannon'")
  a("* Evenness: Shannon_Entropy/log(Richness)")
  a("")
  
  
  # sample_data_table = fread(sample_data_path, data.table = F, select = c(sample_id_column, sample_folder_column))
  if(sample_id_column %ni% names(sample_data_table)){
    error(paste0(sample_id_column, " was not found in sample_data_table."))
  }
  
  if(sample_folder_column %ni% names(sample_data_table)){
    error(paste0(sample_folder_column, " was not found in sample_data_table."))
  }
  
  sample_lut = sample_data_table[[sample_id_column]]
  names(sample_lut) = sample_data_table[[sample_folder_column]]
  
  
  column_output_order = c("Chain", "Count", "aaCDR3", "ntCDR3", "Length", "All_V", "All_J", "Top_V", "Top_J")
  
  # a("Get the files in a common format, not really needed since everything is done with MiXCR, ", 
  #   "but this makes it match up with the diversity extraction step.")
  my_mclapply_output = mclapply(input_file_paths, function(input_file_path){
    output_list = list()
    
    folder_name = basename(dirname(input_file_path))
    sample_name = sample_lut[folder_name]
    
    if(!is.na(sample_name)){
      my_path = input_file_paths[grep(paste0("/", folder_name, "/"), input_file_paths, fixed = T)]
      
      if(length(readLines(my_path)) > 0){
        nt_dt = fread(my_path, data.table = F, 
                      select = c("cloneCount", "aaSeqCDR3", "nSeqCDR3","bestVHit", "bestJHit", 
                                 "allVHitsWithScore", "allJHitsWithScore"))
        nt_dt$Chain = substring(nt_dt$bestVHit , 1, 3)
        names(nt_dt) = c("Count", "aaCDR3","ntCDR3", "Top_V", "Top_J", "All_V", "All_J", "Chain")
        
        nt_dt %<>% run_common_clone_steps()
        
        nt_dt$Top_V = format_all_regions(nt_dt$Top_V)
        nt_dt$Top_J = format_all_regions(nt_dt$Top_J)
        nt_dt$All_V = format_all_regions(nt_dt$All_V, my_split = ",")
        nt_dt$All_J = format_all_regions(nt_dt$All_J, my_split = ",")
        
        nt_dt = nt_dt[ , column_output_order[column_output_order %in% names(nt_dt)]]
        
        
        # now do diversity by aa
        # tapply counts
        aa_counts = tapply(nt_dt$Count, nt_dt$aaCDR3, sum)
        aa_dt = nt_dt[,c("Count", "Chain", "aaCDR3")]
        aa_dt = aa_dt[!duplicated(aa_dt$aaCDR3), ]
        aa_dt$Count = aa_counts[aa_dt$aaCDR3]
        aa_dt = aa_dt[order(aa_dt$Count, decreasing = T), ]
        
        # now do diversity by aa
        # tapply counts
        v_counts = tapply(nt_dt$Count, nt_dt$Top_V, sum)
        v_dt = nt_dt[,c("Count", "Chain", "Top_V")]
        v_dt = v_dt[!duplicated(v_dt$Top_V), ]
        v_dt$Count = v_counts[v_dt$Top_V]
        v_dt = v_dt[order(v_dt$Count, decreasing = T), ]
        
        j_counts = tapply(nt_dt$Count, nt_dt$Top_J, sum)
        j_dt = nt_dt[,c("Count", "Chain", "Top_J")]
        j_dt = j_dt[!duplicated(j_dt$Top_J), ]
        j_dt$Count = j_counts[j_dt$Top_J]
        j_dt = j_dt[order(j_dt$Count, decreasing = T), ]
        
        total_abundance = sum(nt_dt$Count)
        
        message(sample_name)
        output_list["Sample_ID"] = sample_name
        
        # for each chain do our own diversity metrics and write out the vdjtool matrices
        for(my_chain in my_chains){
          chain_df = nt_dt[nt_dt$Chain == my_chain, ]
          my_abundance = sum(chain_df$Count)
          output_list[paste0(my_chain,"_Abundance")] = my_abundance
          output_list[paste0(my_chain,"_Fraction")] = my_abundance/total_abundance
          
          # https://www.biorxiv.org/content/biorxiv/early/2019/06/14/665612.full.pdf
          # TCR convergence was calculated as the aggregate frequency of clones sharing a variable gene 
          # (excluding allele information) and CDR3AA sequence with at least one other identified clone.
          my_richness = nrow(chain_df)
          if(my_richness > 1){
            repeat_aacdr3_counts = summary(factor(chain_df$aaCDR3),  maxsum = Inf)
            output_list[paste0(my_chain,"_Convergence")] = sum(repeat_aacdr3_counts[repeat_aacdr3_counts >1])/my_richness
          } else {
            output_list[paste0(my_chain,"_Convergence")] = NA
          }
          
        }
        
        list_dts = list(nt_dt, aa_dt, v_dt, j_dt)
        names(list_dts) = c("NT", "AA", "VRegion", "JRegion")
        # do nt diversity
        # make individual chain subdat
        # chain_file_name = paste0(my_sample, "_", my_chain, ".tsv")
        for(dt_index in 1:length(list_dts)){
          my_dt = list_dts[[dt_index]]
          my_prefix = names(list_dts)[dt_index]
          
          for(my_chain in my_chains){
            chain_df = my_dt[my_dt$Chain == my_chain, ]
            my_counts = chain_df$Count
            my_abundance = sum(my_counts, na.rm = T)
            chain_richness = sum(my_counts > 0)
            chain_prefix = paste0(my_prefix, "_", my_chain)
            output_list[paste0(chain_prefix,"_Richness")] = chain_richness
            
            if(chain_richness > 1){
              output_list[paste0(chain_prefix,"_d25")] = binfotron::dXX_index(my_counts, 0.25)
              output_list[paste0(chain_prefix,"_d50")] = binfotron::dXX_index(my_counts, 0.50)
              output_list[paste0(chain_prefix,"_d75")] = binfotron::dXX_index(my_counts, 0.75)
              output_list[paste0(chain_prefix,"_Inv_Simpson")] = binfotron::inv_simpson(my_counts)
              output_list[paste0(chain_prefix,"_Chao1")] = binfotron::chao1(my_counts)
              output_list[paste0(chain_prefix,"_Evenness")] = binfotron::evenness(my_counts)
              output_list[paste0(chain_prefix,"_Shannon_Entropy")] = binfotron::shannon_entropy(my_counts)
            } 
          }
        }
      } else {
        # file is there but is empty
        # need to add 0 abundance for all chains and 0 richness for all chians (nt & aa)
        for(my_chain in my_chains){
          output_list[paste0(my_chain,"_Abundance")] = 0
          output_list[paste0("NT_", my_chain,"_Richness")] = 0
          output_list[paste0("AA_", my_chain,"_Richness")] = 0
          output_list[paste0("VRegion_", my_chain,"_Richness")] = 0
          output_list[paste0("JRegion_", my_chain,"_Richness")] = 0
        }
      }
    }
    return(output_list)
  }, mc.cores = 1)
  
  dat = rbindlist(my_mclapply_output, use.names = T, fill = T)
  # now do gene level counts
  
  fwrite(dat, file.path(output_dir, "diversity_data.tsv"), sep = "\t")
}