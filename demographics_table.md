Demographics
================
Ashlyn Johnson

``` r
library(tidyverse)
library(flextable)
library(readxl)
library(gt)
library(gtsummary)
library(officer)
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] officer_0.3.19  gtsummary_1.4.2 gt_0.3.1        readxl_1.3.1   
    ##  [5] flextable_0.6.7 forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7    
    ##  [9] purrr_0.3.4     readr_2.0.1     tidyr_1.1.3     tibble_3.1.4   
    ## [13] ggplot2_3.3.5   tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7          lattice_0.20-44     lubridate_1.7.10   
    ##  [4] assertthat_0.2.1    digest_0.6.27       utf8_1.2.2         
    ##  [7] R6_2.5.1            cellranger_1.1.0    backports_1.2.1    
    ## [10] reprex_2.0.1        evaluate_0.14       httr_1.4.2         
    ## [13] pillar_1.6.2        gdtools_0.2.3       rlang_0.4.11       
    ## [16] uuid_0.1-4          rstudioapi_0.13     data.table_1.14.0  
    ## [19] Matrix_1.3-4        rmarkdown_2.10      splines_4.1.1      
    ## [22] munsell_0.5.0       broom_0.7.9         compiler_4.1.1     
    ## [25] modelr_0.1.8        xfun_0.25           pkgconfig_2.0.3    
    ## [28] systemfonts_1.0.2   base64enc_0.1-3     htmltools_0.5.2    
    ## [31] tidyselect_1.1.1    fansi_0.5.0         crayon_1.4.1       
    ## [34] tzdb_0.1.2          dbplyr_2.1.1        withr_2.4.2        
    ## [37] grid_4.1.1          jsonlite_1.7.2      gtable_0.3.0       
    ## [40] lifecycle_1.0.0     DBI_1.1.1           magrittr_2.0.1     
    ## [43] scales_1.1.1        zip_2.2.0           cli_3.0.1          
    ## [46] stringi_1.7.4       broom.helpers_1.3.0 fs_1.5.0           
    ## [49] xml2_1.3.2          ellipsis_0.3.2      generics_0.1.0     
    ## [52] vctrs_0.3.8         tools_4.1.1         glue_1.4.2         
    ## [55] hms_1.1.0           survival_3.2-11     fastmap_1.1.0      
    ## [58] yaml_2.2.1          colorspace_2.0-2    rvest_1.0.1        
    ## [61] knitr_1.33          haven_2.4.3

### Reading in Data

``` r
samples <-read_csv("raw_data/attributes-file.csv") %>%
  mutate(
    Race = str_to_upper(Race),
    Sex = str_to_upper(Sex),
    Disease = factor(Disease, levels = c("Control", "AD", "FTD"))
  ) %>%
  select(-c(`Sample Name`, `Sample Replacement`)) %>% 
  relocate(c(Disease, `Sample ID`, Age, PMI, APOE, Sex, Race))
```

    ## Rows: 24 Columns: 9

    ## -- Column specification --------------------------------------------------------
    ## Delimiter: ","
    ## chr (7): Sample Name, Sample Replacement, APOE, Race, Sex, Disease, Sample ID
    ## dbl (2): Age, PMI

    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

### Summary Table

``` r
# telling gt summary to calculate mean and sd
theme_gtsummary_mean_sd(set_theme = TRUE)

# creating the summary table
summary_demographics_table <- samples %>%
  select(-`Sample ID`) %>%
  tbl_summary(by = Disease) %>%
  add_p() %>%
  as_flex_table() %>%
  bold(part = "header") %>%
  bold(j = 1)

# save as image
save_as_image(summary_demographics_table, "stats_results/summary_demographics_table.png")
```

    ## [1] "C:/Users/ashly/Documents/HalesLab_AGJ_Local/JohnsonEtAl_GlialProfilingManuscript/stats_results/summary_demographics_table.png"

``` r
save_as_docx(summary_demographics_table, path = "stats_results/summary_demographics_table.docx")
```

![**Summary Demographics
Table**](stats_results/summary_demographics_table.png)
