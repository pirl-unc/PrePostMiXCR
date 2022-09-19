Functions for writing commands and post processing the MiXCR pipeline.

Given a sample data matrix that indicates the sample name, sample folder name and input paths
this script can write out pipeline commands to run MiXCR on the cluster and post process the
data.

## Example commands
In R:
``` r 
diversity_dir = file.path(POST_PROCESSING_DIR, "mixcr")
dir.create(mixcr_dir, showWarnings = F)
output_path = file.path(diversity_dir, "diversity_data.tsv" )


input_file_paths = system(paste0("ls ", RAW_DATA_DIR, "/mixcr/sample_output/*/*_clones.txt"), 
                          intern = TRUE)

my_chains = c( "IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")

if(!("PrePostMiXCR" %in% rownames(installed.packages())))
  devtools::install_github("Benjamin-Vincent-Lab/PrePostMiXCR")
  
PrePostMiXCR::post_process_mixcr(
  input_file_paths = input_file_paths,
  my_chains = my_chains,
  output_path = output_path,
  sample_data_table = fread(sample_data_path),
  thread_num = 1, # works pretty slow with mutliple cores.  probably memory limited.
  sample_folder_column = "run_accession"
)
```


## Assembling this package
In R:
``` r
housekeeping::assemble_package(package_name = "PrePostMiXCR", my_version = "0.0-29",
  my_dir = "/datastore/alldata/shiny-server/rstudio-common/dbortone/packages/PrePostMiXCR")
```

## Push changes
In bash:
``` bash
cd /datastore/alldata/shiny-server/rstudio-common/dbortone/packages/PrePostMiXCR
my_comment="git issues."
git commit -am "$my_comment"; git push origin master
git tag -a 0.0-29 -m "$my_comment"; git push -u origin --tags
```

## Install
Restart R
In R (local library, packrat library):
``` r
devtools::install_github("Benjamin-Vincent-Lab/PrePostMiXCR")
```

Or for a specific version:
``` r
devtools::install_github("Benjamin-Vincent-Lab/PrePostMiXCR", ref = "0.0-27")
```

## Previous location
https://bitbucket.org/unc_lineberger/prepostmixcr/src/master/
