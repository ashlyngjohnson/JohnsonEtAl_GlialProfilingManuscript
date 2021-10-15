Cell Type Score Analysis
================
Ashlyn Johnson

``` r
library(tidyverse)
library(rstatix)
library(ggsci)
library(knitr)
library(officer)
library(flextable)
library(qwraps2)
library(janitor)
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
    ##  [1] janitor_2.1.0   qwraps2_0.5.2   flextable_0.6.7 officer_0.3.19 
    ##  [5] knitr_1.33      ggsci_2.9       rstatix_0.7.0   forcats_0.5.1  
    ##  [9] stringr_1.4.0   dplyr_1.0.7     purrr_0.3.4     readr_2.0.1    
    ## [13] tidyr_1.1.3     tibble_3.1.4    ggplot2_3.3.5   tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.2        jsonlite_1.7.2    carData_3.0-4     modelr_0.1.8     
    ##  [5] assertthat_0.2.1  cellranger_1.1.0  yaml_2.2.1        gdtools_0.2.3    
    ##  [9] pillar_1.6.2      backports_1.2.1   glue_1.4.2        uuid_0.1-4       
    ## [13] digest_0.6.27     rvest_1.0.1       snakecase_0.11.0  colorspace_2.0-2 
    ## [17] htmltools_0.5.2   pkgconfig_2.0.3   broom_0.7.9       haven_2.4.3      
    ## [21] scales_1.1.1      openxlsx_4.2.4    rio_0.5.27        tzdb_0.1.2       
    ## [25] generics_0.1.0    car_3.0-11        ellipsis_0.3.2    withr_2.4.2      
    ## [29] cli_3.0.1         magrittr_2.0.1    crayon_1.4.1      readxl_1.3.1     
    ## [33] evaluate_0.14     fs_1.5.0          fansi_0.5.0       xml2_1.3.2       
    ## [37] foreign_0.8-81    tools_4.1.1       data.table_1.14.0 hms_1.1.0        
    ## [41] lifecycle_1.0.0   munsell_0.5.0     reprex_2.0.1      zip_2.2.0        
    ## [45] compiler_4.1.1    systemfonts_1.0.2 rlang_0.4.11      grid_4.1.1       
    ## [49] rstudioapi_0.13   base64enc_0.1-3   rmarkdown_2.10    gtable_0.3.0     
    ## [53] abind_1.4-5       DBI_1.1.1         curl_4.3.2        R6_2.5.1         
    ## [57] lubridate_1.7.10  fastmap_1.1.0     utf8_1.2.2        stringi_1.7.4    
    ## [61] Rcpp_1.0.7        vctrs_0.3.8       dbplyr_2.1.1      tidyselect_1.1.1 
    ## [65] xfun_0.25

### Read in Data

``` r
celltype_scores <- read_csv("raw_data/cell_type_scores-raw.csv")
```

    ## New names:
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * `` -> ...6
    ## * ...

    ## Rows: 25 Columns: 12

    ## -- Column specification --------------------------------------------------------
    ## Delimiter: ","
    ## chr (12): cell type scores - log2, ...2, ...3, ...4, ...5, ...6, ...7, ...8,...

    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

Based on the advanced analysis

### Initial Munge

``` r
celltype_colnames <- celltype_scores %>%
  slice(1) %>%
  unlist(use.names = FALSE) %>%
  replace_na("sample")

colnames(celltype_scores) <- celltype_colnames

celltype_scores <- celltype_scores %>%
  slice(-1) %>%
  separate("sample", c("Disease", "Case"), sep = "_") %>%
  mutate(Disease = factor(Disease, levels = c("Control", "AD", "FTD"))) %>% 
  mutate_at(c(3:13), as.double)
```

### ANOVA

``` r
anova_results <- sapply(celltype_scores[, -c(1:2)], function(x) {
  lm(x ~ Disease, data = celltype_scores) %>%
    anova() %>%
    unlist()
})

adjusted_anova_results <- anova_results %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Cell Type") %>%
  select("Cell Type", Fval = "F value1", pval = "Pr(>F)1") %>%
  mutate(p.adjval = p.adjust(pval, method = "BY")) 

significant_celltypes <- adjusted_anova_results %>% 
  filter(p.adjval <= .05) %>% 
  pull("Cell Type")

kable(adjusted_anova_results, caption = "All One Way ANOVA Results with Correction for Multiple Tests")
```

| Cell Type              |      Fval |      pval |  p.adjval |
|:-----------------------|----------:|----------:|----------:|
| Astrocytes             | 11.636416 | 0.0003971 | 0.0131899 |
| Neutrophils            |  6.024607 | 0.0085528 | 0.0757888 |
| Oligodendrocytes       |  7.559545 | 0.0033656 | 0.0559008 |
| B-cells                |  2.358716 | 0.1191003 | 0.4395946 |
| Cytotoxic cells        |  3.388818 | 0.0530276 | 0.2516436 |
| Mast cells             |  2.153454 | 0.1410237 | 0.4684617 |
| Endothelial Cells      |  3.149874 | 0.0636264 | 0.2641977 |
| DC                     |  3.771097 | 0.0398736 | 0.2207578 |
| Macrophages\_Microglia |  5.429078 | 0.0125744 | 0.0835407 |
| Neuron                 |  1.351354 | 0.2804942 | 0.8470581 |
| CD45                   |  5.922833 | 0.0091261 | 0.0757888 |

All One Way ANOVA Results with Correction for Multiple Tests

### Pairwise T-Tests

``` r
celltype_scores_tidy <- celltype_scores %>%
  select(one_of(c("Disease", "Case", significant_celltypes))) %>%
  pivot_longer(3,
    names_to = "celltype",
    values_to = "celltype_score"
  )

posthoc_t_test_results <- celltype_scores_tidy %>%
  group_by(celltype) %>%
  pairwise_t_test(data = ., celltype_score ~ Disease, pool.sd = TRUE, p.adjust = "none") %>%
  select(1:8) %>%
  mutate(p.adj = p.adjust(p, method = "BY"))

sig_posthoc_t_test_results <- posthoc_t_test_results %>% 
  filter(p.adj <=.05)

kable(posthoc_t_test_results, caption = "Pairwise T-Test Results for Pathways with Significant ANOVA P-Vals")
```

| celltype   | .y.             | group1  | group2 |  n1 |  n2 |        p | p.signif |     p.adj |
|:-----------|:----------------|:--------|:-------|----:|----:|---------:|:---------|----------:|
| Astrocytes | celltype\_score | Control | AD     |   8 |   8 | 5.90e-02 | ns       | 0.1081667 |
| Astrocytes | celltype\_score | Control | FTD    |   8 |   8 | 9.59e-05 | \*\*\*\* | 0.0005275 |
| Astrocytes | celltype\_score | AD      | FTD    |   8 |   8 | 1.06e-02 | \*       | 0.0291500 |

Pairwise T-Test Results for Pathways with Significant ANOVA P-Vals

### Creating Results File

``` r
write_csv(adjusted_anova_results, "stats_results/celltype_analysis_anova.csv")
write_csv(posthoc_t_test_results, "stats_results/celltype_analysis_posthoc.csv")
```

### Plots with P-vals

Here, I want draw plots for the pathways that include ANOVAs p values.

``` r
# need a vector with all the cell types
celltypes <- colnames(celltype_scores)[3:13]


# loop to generate plots
for (i in seq_along(celltypes)) {
  plot <- celltype_scores %>%
    ggplot(aes(x = Disease, y = .data[[celltypes[i]]])) +
    geom_boxplot(aes(fill = Disease)) +
    scale_fill_manual(values = c("turquoise2", "magenta4", "yellow3" )) +
    geom_point(size = 3) +
    ggtitle(celltypes[i]) +
    ylab("Cell Type Score") +
    labs(subtitle = paste0("One Way ANOVA, adjusted p-value = ", adjusted_anova_results %>%
      filter(.data[["Cell Type"]] == celltypes[i]) %>%
      select(p.adjval) %>%
      round(digits = 5))) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13)
    )

  ggsave(paste0("figures/celltype_scores/", str_replace_all(celltypes[i], " ", ""), ".png"))
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

### Nice Table

``` r
celltype_scores_summary_table <- celltype_scores %>%
  select(-Case) %>%
  group_by(Disease) %>%
  summarize(across(where(is.double), ~ mean_sd(.x, denote_sd = "paren"))) %>%
  t() %>%
  as.data.frame() %>%
  row_to_names(row_number = 1) %>%
  rownames_to_column() %>%
  rename("Cell Type" = rowname) %>%
  left_join(adjusted_anova_results, by = "Cell Type") %>% # joining with anova results
  rename("F-value" = Fval,
         "P-value" = pval,
  ) %>%
  mutate(
    "F-value" = round(.data[["F-value"]], digits = 3),
    "P-value" = round(.data[["P-value"]], digits = 5),
    p.adjval = round(p.adjval, digits = 5)
  ) %>%
flextable() %>% # piping into flextable
  add_header_row( # creating another header row to label mean(sd) and anova results
    values = c(
      "",
      "Mean (SD)",
      "Mean (SD)",
      "Mean (SD)",
      "ANOVA Results",
      "ANOVA Results",
      "ANOVA Results"
    )
  ) %>%
  merge_h(part = "header") %>%
  align(align = "center", part = "header") %>%
  align(align = "Center", part = "header") %>%
  vline( # creating vertical lines
    j = c("Cell Type", "FTD"),
    border = fp_border(color = "gray"),
    part = "body"
  ) %>%
  bold(part = "header") %>%
  fontsize(i = 1, size = 14, part = "header") %>%
  fontsize(i = 2, size = 13, part = "header") %>%
  fontsize(size = 12, part = "body") %>%
  # bg(i = ~ p.adjval <= .05, # highlighting the significant results 
  #    bg = "lightblue",
  #    part = "body") %>%
  set_header_labels(p.adjval = "Adj P-value")

# saving as an image

save_as_image(celltype_scores_summary_table, "stats_results/celltype_scores_table.png")
save_as_pptx(celltype_scores_summary_table, path = "stats_results/celltype_scores_table.pptx")
```

![**Cell Type Scores**](stats_results/celltype_scores_table.png)
