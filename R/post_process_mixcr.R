
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' post_process_star_salmon
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title post_process_star_salmon 
#' 
#' @description 
#' Joins individual sample files into one tsv file.  Also converts ucsc names to hgnc_symbol|entrez_ids
#' and outputs the gene level counts.
#' 
#' @param input_file_paths Character vector of paths to the pipeline output data. 
#' @param output_dir Path to the output folder.
#' @param sample_data_path Path to the sample data which should contatin the sample_id_column and sample_folder_column
#' @param gene_biotypes The type of biomaRt gene_biotypes that should be output to the gene level output.
#' @param sample_folder_column The name of the column that has sample folder names
#' @param sample_id_column The name of the column that has sample ids
#' @param thread_num Integer number of threads to run mclapply statements
#' 
#' @return A path to the rds file.
#' 
#' @family mart
#' 
#' @export
post_process_mixcr = function(
  # this_script_path = get_script_dir_path(include_file_name = T)
  # init_path = find_file_along_path(this_script_path, "_init.R")
  # source(init_path)
  # base_dir = dirname(init_path)
  
  input_file_paths,# = system(paste0("ls ", RAW_DATA_DIR, "/pipeline_output/star_salmon/*/*_quant.sf"), intern = TRUE)
  output_dir,# = file.path(base_dir, "post_processing", "star_salmon")
  sample_data_path,# = file.path(base_dir, "sample_data", "sample_data.tsv")
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
  a("NT metrics are by nucleotide sequence diversity")
  a("AA metrics are by amino acid sequence diversity")
  a("VRegion metrics are by V region (top hit) diversity")
  a("JRegion metrics are by J region (top hit) diversity")
  a("Chain abundance and fraction are the same across NT, AA, VRegion, and JRegion")
  a("")
  
  sample_dat = fread(sample_data_path, data.table = F, select = c(sample_id_column, sample_folder_column))
  sample_lut = sample_dat[[sample_id_column]]
  names(sample_lut) = sample_dat[[sample_folder_column]]
  
  
  column_output_order = c("Chain", "Count", "aaCDR3", "ntCDR3", "Length", "All_V", "All_J", "Top_V", "Top_J")
  
  # a("Get the files in a common format, not really needed since everything is done with MiXCR, ", 
  #   "but this makes it match up with the diversity extraction step.")
  my_mclapply_output = mclapply(input_file_paths, function(input_file_path){
    
    folder_name = basename(dirname(input_file_path))
    sample_name = sample_lut[folder_name]
    
    my_path = input_file_paths[grep(paste0("/", folder_name, "/"), input_file_paths, fixed = T)]
    
    nt_dt = fread(my_path, data.table = F, 
                  select = c("cloneCount", "Chains", "aaSeqCDR3", "nSeqCDR3","bestVHit", "bestJHit", 
                             "allVHitsWithScore", "allJHitsWithScore"))
    names(nt_dt) = c("Count", "Chain", "aaCDR3","ntCDR3", "Top_V", "Top_J", "All_V", "All_J")
    
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
     
    output_list = list()
    message(sample_name)
    output_list["Sample_ID"] = sample_name
    
    # for each chain do our own diversity metrics and write out the vdjtool matrices
    for(my_chain in my_chains){
      chain_df = nt_dt[nt_dt$Chain == my_chain, ]
      my_abundance = sum(chain_df$Count)
      output_list[paste0(my_chain,"_Abundance")] = my_abundance
      output_list[paste0(my_chain,"_Fraction")] = my_abundance/total_abundance
      
      #' https://www.biorxiv.org/content/biorxiv/early/2019/06/14/665612.full.pdf
      #' TCR convergence was calculated as the aggregate frequency of clones sharing a variable gene 
      #' (excluding allele information) and CDR3AA sequence with at least one other identified clone.
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
        } else {
          output_list[paste0(chain_prefix,"_d25")] = NA
          output_list[paste0(chain_prefix,"_d50")] = NA
          output_list[paste0(chain_prefix,"_d75")] = NA
          output_list[paste0(chain_prefix,"_Inv_Simpson")] = NA
          output_list[paste0(chain_prefix,"_Chao1")] = NA
          output_list[paste0(chain_prefix,"_Evenness")] = NA
          output_list[paste0(chain_prefix,"_Shannon_Entropy")] = NA
        }
      }
    }
    return(output_list)
  }, mc.cores = 1)
  
  dat = rbindlist(my_mclapply_output, use.names = T, fill = T)
  # now do gene level counts
  
  fwrite(dat, file.path(output_dir, "diversity_data.tsv"), sep = "\t")
}