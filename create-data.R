


## COMMON SOURCE FOR ALL SHINY KIDNEY APPS
   ## Based on R_ccpRCC_2019_SOURCE.R # 2020-05-12


## ==============================
#        SOURCE CODE
## ==============================

require(dlfoo2)
require(dlfoo2data)
require(Biobase)


require(reshape2)
require(dplyr)
options(dplyr.width = Inf)
require(tibble)
require(tidyverse)
require(tidyr)

setwd("~/OneDrive - Lund University/SHINY_KIDNEY/")


## :::::::::::::
##    LOAD DATA
## ::::::::::::::
## load data stored in dlfoo2 package
   #data("tcga_surv_tab", package="dlfoo2data"
   data("shape_subtypes", package="dlfoo2")
    data("pan_gdc", package="dlfoo2data")




# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#      PDATA SOURCE
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

pdata_source <- as.data.frame(dlfoo2::pdata_panCanFull ) %>% dplyr::filter(type %in% "BRCA")
pdata_source <- pdata_source %>% droplevels()
   # %>% mutate(tax_simp3 = as.character(tax_simp2)) %>%  mutate(tax_simp3 = if_else(is.na(tax_simp3) & grepl("-NT", sample_type2), "kidney_n", tax_simp3)) %>% mutate(tax_simp3 = factor(tax_simp3, levels=c("kidney_n", levels(tax_simp2))))
   # tax_simp3 = all normals annotated
u <- grepl("-NT", pdata_source$sample_type2)
pdata_source$purity[u] <- NA
pdata_source$ploidy[u] <- NA
pdata_source$Cancer.DNA.fraction[u] <- NA
pdata_source$Subclonal.genome.fraction[u] <- NA


   pdata_temp <- pdata_source %>% dplyr::filter(platform =="IlluminaHiSeq_RNASeqV2") %>% dplyr::filter(sample_id %in% sampleNames(pan_gdc)) %>%
      dplyr::filter(!duplicated(sample_id))
   str(pdata_temp)
   my_samples <- pdata_temp$sample_id

   pdata_jonas <- read.delim("Data/BRCA_clinicalMatrix.txt", as.is=T)
   str(pdata_jonas)
   table(duplicated(pdata_jonas$sampleID))
   # sampleID "CGA-3C-AAAU-01""

   my_samples <- intersect(pdata_temp$sample_barcode, pdata_jonas$sampleID)

   pdata_brca <- pdata_jonas %>% dplyr::filter(sampleID %in% my_samples)
   pdata_temp <- pdata_temp %>% dplyr::rename(sampleID = sample_barcode)
   pdata_brca <- pdata_brca %>%
      dplyr::left_join(pdata_temp %>% dplyr::select(sample_id, sampleID))

      brca_nt <- pan_gdc[,pdata_brca$sample_id]

   identical(sampleNames(brca_nt), pdata_brca$sample_id)


   saveRDS(fData(brca_nt), "Data/jonas_brca_gene_table.rds")
   saveRDS(exprs(brca_nt), "Data/jonas_brca_gene_expr_matrixe.rds")
   saveRDS(pdata_brca, "Data/jonas_brca_sample_sheet.rds")

   write.table(data.frame(annot=colnames(pdata_brca)), sep="\t", quote=F, row.names=F, file="Data/jonas_pdata_columns.txt")

#
#
# # get pancan tumor color codes
# pancan_colors <- read.delim("Data/PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv", as.is=T)
# u <- match(pancan_colors$Study.Abbreviation, names(color_subtypes))
# color_subtypes[u] <- pancan_colors$Hex.Colors
#
# temp_vec <- c("#006699","goldenrod2")
# names(temp_vec) <- c("KIRC-TAP", "KIRP-TAP")
# color_subtypes <- c(color_subtypes, temp_vec)
#
# myfile <- "Data/retax_2020_colorTable.txt"
# retax_colors <- read.delim(myfile, as.is=T)$color
# names(retax_colors) <- read.delim(myfile, as.is=T)$name
# color_subtypes <- c(retax_colors, color_subtypes)
# color_subtypes["Rest"] <- "gray"
#
# match(levels(pdata_source$tax), names(color_subtypes))
#
# retax_shape <- read.delim("Data/retax_2020_colorTable.txt", as.is=T)$shape
# names(retax_shape) <- read.delim("Data/retax_2020_colorTable.txt", as.is=T)$name
# shape_subtypes <- c(retax_shape, shape_subtypes)
#
#
# ## color_annotations
# color_subtypes <- c(color_subtypes, color_annotations)
# # uu <- unique(unlist(apply((pdata_bin %>% select(-sample_id)),2,unique)))
# # u<-match(uu, c(names(color_annotations), names(color_subtypes)))
# # u
# # uu[is.na(u)]





## TRANSCRIPTION FACTORS
## ---------------------

tf_df <- read.delim("Data/Human Transcription Factors Cell 2018.txt", as.is=T)
head(tf_df)
tf_df <- tf_df %>% dplyr::rename(ENSG=ID, SYMBOL=Name) %>% dplyr::filter(is.TF=="Yes")
rownames(tf_df) <- tf_df$ENSG
my_ensg <- intersect(rownames(tf_df), featureNames(rcc_nt))
tf_df <- tf_df[my_ensg,]
tf_mat <- centroid_list$tax_simp2[tf_df$ENSG, levels(pdata_rccGex$tax_simp3)]
tf_df$lfc.taxonomy.maxvar <- apply(tf_mat, 1, function(x) round(max(x)-min(x),2))
head(tf_df)


## each group vs next highest group
my_groups <- c(levels(pdata_source$tax_simp3))
cent_mat <- centroid_list$tax_simp2[rownames(tf_mat), my_groups]
group_diff <- apply(cent_mat, 1, function(x){
   yy <- vector(length=8L, mode="numeric")
   for(i in 1:length(x)){
      y <- x[-i]
      #yy[i] <- round(x[i]-median(y),2)
      yy[i] <- round(x[i]-max(y),2)
   }
   return(yy)
})
group_diff <- t(group_diff)
colnames(group_diff) <- paste0(my_groups,".vs.hi")
tf_df <- cbind(tf_df, group_diff)
tf_df <- tf_df %>% select(ENSG, SYMBOL, lfc.taxonomy.maxvar, !!!colnames(group_diff), everything()) %>% arrange(-lfc.taxonomy.maxvar)
head(tf_df)
rm(my_groups, cent_mat, tf_mat, my_ensg, group_diff)





### LIMMA
