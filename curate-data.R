



## Limma results table
## -----------------
 # - table with contrast, p.value
 # Log fiold difference from centroid_list


library(shiny)
setwd("~/OneDrive - Lund University/SHINY_KIDNEY/")
source("R_SOURCE_SHINY_KIDNEY.R")

source("R_SOURCE_SHINY_KIDNEY_FUNCTIONS.R")


#
# pdata = pdata_source
# gene_df <- fData(rcc_nt) %>% dplyr::select(ENSG, SYMBOL,ENTREZID, chr)
#
# annotation_bin_choices <- colnames(pdata_bin %>% dplyr::select(-sample_id))
# default.annotation_bin <- "tax_simp"
# annotation_num_choices <- colnames(pdata_num %>% dplyr::select(-sample_id))
# default.annotation_num<- "HNF_score"
#
#

# ## For Limma
# ## ---------
#    my_fit <- readRDS(file="Data/Limma_taxSimp2_fit2.rds")
#    contrast_tokens <- colnames(my_fit$coefficients)
#    u <- unlist(strsplit(contrast_tokens, split = " - "))
#    contrast_1_tokens <- unique(u[seq_along(u) %% 2 == 1])
#    contrast_ref_tokens  <- unique(u[seq_along(u) %% 2 == 0])
#
#    t_df <- reshape2::melt(my_fit$t, variable.name="ENSG", value.name="t.value") %>%
#       dplyr::rename(contrast = "Contrasts", ENSG="Var1") %>%
#       mutate(sample_id = paste(contrast, ENSG))
#    str(t_df)
#    x_df <- reshape2::melt(my_fit$p.value, variable.name="ENSG", value.name="p.value") %>%
#          dplyr::rename(contrast = "Contrasts", ENSG="Var1") %>%
#          mutate(sample_id = paste(contrast, ENSG))
#    x_df <- dplyr::left_join(t_df, x_df %>% dplyr::select(-ENSG, -contrast), by="sample_id")
#    str(x_df)
#
#    x_df$SYMBOL <- my_fit$genes$SYMBOL[match(x_df$ENSG, my_fit$genes$ENSG)]
#    x_df$c1 <- unlist(sapply(as.character(x_df$contrast), function(x) strsplit(x, split = " - ")[[1]][1]))
#    x_df$c2 <- unlist(sapply(as.character(x_df$contrast), function(x) strsplit(x, split = " - ")[[1]][2]))
#    x_df$lfc <- unlist(lapply(as.character(unique(x_df$contrast)), function(x){
#       centroid_list$tax_simp2[, strsplit(x, split = " - ")[[1]][1]] -  centroid_list$tax_simp2[,strsplit(x, split = " - ")[[1]][2]]
#          }))
#    message(" .. calculating BH adjusted p.vals for all contrasts")
#    x_df$gene_padj_BH_call = as.character(unlist(lapply(as.character(unique(x_df$contrast)), function(x) decideTests(my_fit,  p.value = 0.05, adjust.method = "BH")[,x])))
#    x_df <- x_df %>%
#       mutate(gene_padj_BH_call = dplyr::recode(as.character(gene_padj_BH_call), "0"="ns", "-1"="Downreg","1"="Upreg"))
#    limma_res_df <- x_df
#    rm(x_df)
#    str(limma_res_df)
#
# # limma_res_df %>% mutate_at(vars(c1,c2,gene_padj_BH_call,SYMBOL), funs(as.factor))
# saveRDS(limma_res_df, "Data/limma_t_p_BH_results_table.rds")





## For Limma -- Mahzaar
## ---------
my_fit <- readRDS(file="Data/limma_rlog_transformed.rds")
contrast_tokens <- colnames(my_fit$coefficients)
u <- unlist(strsplit(contrast_tokens, split = " - "))
contrast_1_tokens <- unique(u[seq_along(u) %% 2 == 1])
contrast_ref_tokens  <- unique(u[seq_along(u) %% 2 == 0])

t_df <- reshape2::melt(my_fit$t, variable.name="ENSG", value.name="t.value") %>%
   dplyr::rename(contrast = "Contrasts", ENSG="Var1") %>%
   mutate(sample_id = paste(contrast, ENSG))
str(t_df)
x_df <- reshape2::melt(my_fit$p.value, variable.name="ENSG", value.name="p.value") %>%
   dplyr::rename(contrast = "Contrasts", ENSG="Var1") %>%
   mutate(sample_id = paste(contrast, ENSG))
x_df <- dplyr::left_join(t_df, x_df %>% dplyr::select(-ENSG, -contrast), by="sample_id")
str(x_df)

x_df$SYMBOL <- my_fit$genes$SYMBOL[match(x_df$ENSG, my_fit$genes$ENSG)]
x_df$c1 <- unlist(sapply(as.character(x_df$contrast), function(x) strsplit(x, split = " - ")[[1]][1]))
x_df$c2 <- unlist(sapply(as.character(x_df$contrast), function(x) strsplit(x, split = " - ")[[1]][2]))
x_df$lfc <- unlist(lapply(as.character(unique(x_df$contrast)), function(x){
   centroid_mazhar[, strsplit(x, split = " - ")[[1]][1]] -  centroid_mazhar[,strsplit(x, split = " - ")[[1]][2]]
}))
message(" .. calculating BH adjusted p.vals for all contrasts")
x_df$gene_padj_BH_call = as.character(unlist(lapply(as.character(unique(x_df$contrast)), function(x) decideTests(my_fit,  p.value = 0.05, adjust.method = "BH")[,x])))
x_df <- x_df %>%
   mutate(gene_padj_BH_call = dplyr::recode(as.character(gene_padj_BH_call), "0"="ns", "-1"="Downreg","1"="Upreg"))
limma_res_df <- x_df
rm(x_df)
str(limma_res_df)

# limma_res_df %>% mutate_at(vars(c1,c2,gene_padj_BH_call,SYMBOL), funs(as.factor))
saveRDS(limma_res_df, "Data/limma_t_p_BH_results_table_Mahzar.rds")
