Functions for writing commands and post processing the MiXCR pipeline.

Given a sample data matrix that indicates the sample name, sample folder name and input paths
this script can write out pipeline commands to run MiXCR on the cluster and post process the
data.

Erase this: Test jira tagging a commit message 5.

## Assembling this package
In R:
``` r
housekeeping::assemble_package(package_name = "PrePostMiXCR", my_version = "0.0-20",
  my_dir = "/datastore/alldata/shiny-server/rstudio-common/dbortone/packages/pre_post_mixcr")
```

## Push changes
In bash:
``` bash
cd /datastore/alldata/shiny-server/rstudio-common/dbortone/packages/pre_post_mixcr
my_comment="Renamed so it won't be confused with the docker image repo."
git commit -am "$my_comment"; git push origin master
git tag -a 0.0-13 -m "$my_comment"; git push -u origin --tags
```

## Install
Restart R
In R (local library, packrat library):
``` r
devtools::install_bitbucket("unc_lineberger/pre-post-mixcr")
```

Or for a specific version:
``` r
devtools::install_bitbucket("unc_lineberger/pre-post-mixcr", ref = "0.0-12")
```