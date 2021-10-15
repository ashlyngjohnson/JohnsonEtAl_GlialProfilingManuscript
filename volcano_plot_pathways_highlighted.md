Volcano Plots with Pathways Highlighted for FTD vs Control
================
Ashlyn Johnson

``` r
library(tidyverse)
library(ggrepel)
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
    ##  [1] ggrepel_0.9.1   forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7    
    ##  [5] purrr_0.3.4     readr_2.0.1     tidyr_1.1.3     tibble_3.1.4   
    ##  [9] ggplot2_3.3.5   tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.1 xfun_0.25        haven_2.4.3      colorspace_2.0-2
    ##  [5] vctrs_0.3.8      generics_0.1.0   htmltools_0.5.2  yaml_2.2.1      
    ##  [9] utf8_1.2.2       rlang_0.4.11     pillar_1.6.2     glue_1.4.2      
    ## [13] withr_2.4.2      DBI_1.1.1        dbplyr_2.1.1     modelr_0.1.8    
    ## [17] readxl_1.3.1     lifecycle_1.0.0  munsell_0.5.0    gtable_0.3.0    
    ## [21] cellranger_1.1.0 rvest_1.0.1      evaluate_0.14    knitr_1.33      
    ## [25] tzdb_0.1.2       fastmap_1.1.0    fansi_0.5.0      broom_0.7.9     
    ## [29] Rcpp_1.0.7       scales_1.1.1     backports_1.2.1  jsonlite_1.7.2  
    ## [33] fs_1.5.0         hms_1.1.0        digest_0.6.27    stringi_1.7.4   
    ## [37] grid_4.1.1       cli_3.0.1        tools_4.1.1      magrittr_2.0.1  
    ## [41] crayon_1.4.1     pkgconfig_2.0.3  ellipsis_0.3.2   xml2_1.3.2      
    ## [45] reprex_2.0.1     lubridate_1.7.10 rstudioapi_0.13  assertthat_0.2.1
    ## [49] rmarkdown_2.10   httr_1.4.2       R6_2.5.1         compiler_4.1.1

### Read in Data

``` r
# FTD diffex info
diffex_FTD <- read_csv("raw_data/DEResults-DiseaseFTD.csv") %>% 
  rename(Gene = 1) %>% 
  mutate(Gene = str_remove_all(Gene, "-mRNA")) %>% 
  mutate(test = BY.p.value < .05)

# importing pathway analysis results to get a list of statistically significant pathways 
sig_pathways <- read_csv("stats_results/pathway_analysis_anova.csv") %>% 
  filter(p.adjval <= .05) %>% 
  pull(pathway) 
```

### Volcano Plots with Pathway Ppecific Genes as a Different Shape

``` r
for (i in seq_along(sig_pathways)) {
  ggplot() +
    geom_point(
      data = diffex_FTD %>%
        filter(str_detect(Gene.sets, sig_pathways[i]) == FALSE),
      aes(
        x = `Log2 fold change`,
        y = -log10(`P-value`),
        color = test
      ),
      size = 2,
      alpha = .5
    ) +
    scale_color_manual(values = c("gray", "yellow3")) +
    geom_point(
      data = diffex_FTD %>%
        filter(str_detect(Gene.sets, sig_pathways[i])),
      aes(
        x = `Log2 fold change`,
        y = -log10(`P-value`),
        color = test
      ),
      size = 3,
      alpha = 1,
      shape = 17
    ) +
    labs(
      y = expression(-log[10]("p-value")),
      x = expression(log[2]("FoldChange")),
      title = paste0(
        "Differential Expression in FTD-tau vs Control:\n",
        sig_pathways[i]
      )
    ) +
    guides(color = guide_legend(reverse = TRUE)) +
    geom_text_repel(
      data = diffex_FTD %>%
        filter(str_detect(Gene.sets, sig_pathways[i]),
               test == TRUE),
      aes(
        x = `Log2 fold change`,
        y = -log10(`P-value`),
        label = Gene
      ),
      size = 4,
      color = "dark red"
    ) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'none')
  
  ggsave(
    paste0(
      "figures/volcanoandpathways/volcano_ftdvsctrl_",
      str_replace_all(sig_pathways[i], " ", ""),
      ".png"
    )
  )
  
}
```

    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image

    ## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
    ## Saving 7 x 5 in image
