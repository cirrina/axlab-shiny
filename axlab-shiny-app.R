# R -e "shiny::runApp('~/scripts/axlab/SHINY_AXLAB/shinyAxlab_app.R', launch.browser = TRUE)"


## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##   SOURCE FILES
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(shiny)
setwd("~/Scripts/axlab/SHINY_AXLAB/")

source("R_Source.R")
# source("R_SOURCE_SHINY_KIDNEY_FUNCTIONS.R")




## expressionSets
## ----------------------------------------------
eset_list <- list(
   rcc_t = rcc_t,
   rcc_nt = rcc_nt,
   pan_gdc = pan_gdc,
   ccle  = ccle,
   mazhar = es_mazhar
)


## Lists of available data sets, fdata, and pdata
## ----------------------------------------------
## This is used to select avtive (reactive) data matrix, fData(eset_list[[input$dataset]]) fdata and pdata

dataset_selection_vec <- names(eset_list)


## a list for the selection of genes - list of gene symobl -> matched ENSG
gene_selection_list <- lapply(eset_list, function(x){
   split(fData(x)[,"ENSG"],fData(x)[,"SYMBOL"])
})


## pdata_list
pdata_list <- lapply(eset_list, function(x){
   pData(x)
})

## List of co;lumns for each data set split in bin/num depending if numeric or not
pdata_column_types_list <- lapply(eset_list, function(x){
   p <- pData(x)
   u <- sapply(colnames(p), function(z) is.numeric(p[,z]))
   return(list(
      bin=names(u)[!u],
      num=names(u)[u]
   ))
})






### sample_id
## -----------------

# sample_id_list <- list()
# for(i in names(pdata_list)){
#    # sample_id_choices <- split(pdata_sample_id[,"sample_id"], pdata_sample_id[,"taxonomy"])
#    sample_id_list[[i]] <- pdata_list[[i]]$sample_id
# }



# LIMMA RESULTS TABLE / not available for all datasets
## --------------------------------------
limmma_res_list <- list(
   rcc_t = readRDS("Data/Limma_taxSimp2_allGenes__t_p_BH_ResultsTable.rds"),
   rcc_nt =readRDS("Data/Limma_taxSimp2_allGenes__t_p_BH_ResultsTable.rds"),
   mazhar = readRDS("Data/limma_t_p_BH_results_table_Mahzar.rds")
)

fdata_limma_list <- list()
contrast_token_list <- list()
contrast_1_token_list <- list()
contrast_ref_token_list <- list()

for(i in names(limmma_res_list)){
   fdata_limma_list[[i]] <- limmma_res_list[[i]][!duplicated(limmma_res_list[[i]]$ENSG),] %>% dplyr::select(ENSG,SYMBOL) #  feature data of all genes in limma data
   contrast_token_list[[i]] <- as.character(unique(limmma_res_list[[i]]$contrast))
   u <- unlist(strsplit(contrast_token_list[[i]], split = " - "))
   contrast_1_token_list[[i]] <- unique(u[seq_along(u) %% 2 == 1])
   contrast_ref_token_list[[i]]  <- unique(u[seq_along(u) %% 2 == 0])
   }






#
# ## ------------------------------------------
# ##  GENE CENTROID MATRIXES (RCC taxonomies)
# ## ------------------------------------------
#
# ## Calcuate gex centroid matrices for all groups - used for log2 diff for e.g. contrasts
# ## changed 20200916 according no new tax203
# ## ---------------
# table(pdata_rccGex$tax)
#
# table(pdata_rccGex$tax, useNA="always")
# table(pdata_rccGex$tax_simp2, useNA="always")
# pdata_rccGex[is.na(pdata_rccGex$tax), ]
#
#
# my_annots <-  c("tax","tax_simp2","tax_simp3","taxonomy_published")
# centroid_list <- lapply(my_annots, function(x){
#    pdata_temp <- pdata_rccGex %>%
#       dplyr::rename(tax_centroid = c(x)) %>%
#       dplyr::filter(!is.na(tax_centroid)) %>%
#       dplyr::filter(!tax_centroid %in% c("MA","MTSCC","potential_swap")) %>%
#       droplevels()
#
#    #colnames(pdata_temp)[2] <- "annot"
#    my.es <- rcc_nt[,pdata_temp$sample_id]
#    pData(my.es) <- pdata_temp
#    table(pdata_temp$tax_centroid)
#    my_mat <- esetGroupMeans(my.es, annot = pdata_temp$tax_centroid)
#    my_mat <- t(apply(my_mat, 1, function(x) x-mean(x)))
# })
# names(centroid_list) <- my_annots
#
#
# my_annots <-  c("tax_simp3")
# centroid_list_notRowCentered <- lapply(my_annots, function(x){
#    pdata_temp <- pdata_rccGex %>%
#       dplyr::rename(tax_centroid = c(x)) %>%
#       dplyr::filter(!is.na(tax_centroid)) %>%
#       dplyr::filter(!tax_centroid %in% c("MA","MTSCC")) %>%
#       droplevels()
#
#    #colnames(pdata_temp)[2] <- "annot"
#    my.es <- rcc_nt[,pdata_temp$sample_id]
#    pData(my.es) <- pdata_temp
#    table(pdata_temp$tax_centroid)
#    my_mat <- esetGroupMeans(my.es, annot = pdata_temp$tax_centroid)
#    #my_mat <- t(apply(my_mat, 1, function(x) x-mean(x)))
# })
# names(centroid_list_notRowCentered) <- my_annots
#



## TF centroids
# tf_mat <- centroid_list$tax_simp2[tf_df$ENSG, ]
# tf_df$lfc.taxonomy.maxvar <- apply(tf_mat, 1, function(x) round(max(x)-min(x),2))
# head(tf_df)
# ## each group vs next highest group
# #my_groups <- c(levels(pdata_source$tax_simp3))
# cent_mat <- centroid_list$tax_simp3[rownames(tf_mat), ]
# group_diff <- apply(cent_mat, 1, function(x){
#    yy <- vector(length=8L, mode="numeric")
#    for(i in 1:length(x)){
#       y <- x[-i]
#       #yy[i] <- round(x[i]-median(y),2)
#       yy[i] <- round(x[i]-max(y),2)
#    }
#    return(yy)
# })
# group_diff <- t(group_diff)
# tf_df <- cbind(tf_df, group_diff)
# tf_df <- tf_df %>% select(ENSG, SYMBOL, lfc.taxonomy.maxvar, !!!colnames(group_diff), everything()) %>% arrange(-lfc.taxonomy.maxvar)
# head(tf_df)
# rm(cent_mat, tf_mat, my_ensg, group_diff)
#
#





## START SHINY
message("## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
message("##   GENE EXPLORER APP WITH TABS ")
message("## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")

rm(ui, server)


ui <- fluidPage(
   titlePanel("SHINY AXLAB"),
      sidebarLayout(
         # Sidebar panel for inputs ----
          sidebarPanel(width=2,
            uiOutput("dataset_selection"),
            uiOutput("gene1_selection"),
            uiOutput("gene2_selection"),
            # uiOutput("sample_id_selection"),
            #  tags$hr(),
            uiOutput("annotation_bin_selection"),
            uiOutput("annotation_num_selection"),
            # # tags$hr(),
            # # uiOutput("contrast_1_selection"),
            # # uiOutput("contrast_2_selection"),
            # tags$hr(),
            tags$hr()
         ),

      # Main panel for displaying outputs ----
       mainPanel(
         # Output: Tabset
         tabsetPanel(type = "tabs",


         ## :::::::::::::::::::::::::::::::::::::
         ##          TAB
         ## :::::::::::::::::::::::::::::::::::::
         ##    Primary Gene plot tab - RCC
         ## :::::::::::::::::::::::::::::::::::::
            ## RCC-centered plots (not PAN tcga)
            ## Violins of 1st and 2nd gene

            tabPanel("Basic Gene Plots - RCC",
               ## Selected gene info
               verbatimTextOutput('gene1_selection_text_t1'),
               fluidPage(htmlOutput("html_link")),
               tags$hr(),
               fluidPage(
                  verbatimTextOutput('')
                  ),
               ## Plot gene 1 vs primary discrete (binned) annotation
               tags$hr(),
               fluidPage(
                  column(8, plotlyOutput('violin_annot_bin_gene1'), height=1200)
                  ),
               # ## ScatterPlot gene 1 vs Primary numeric annotation
               # ## ScatterPlot gene 1 vs gene 2
                fluidPage(
                    column(6, plotlyOutput('scatter_gene1_gene2'), height=1200),
                    column(6, plotlyOutput('scatter_annot_num_gene1'), height=1200)
                   )
               ),


         ## :::::::::::::::::::::::::::::::::::::
         ##              TAB
         ## :::::::::::::::::::::::::::::::::::::
         ##       GENE SIGNATURES KIDNEY
         ## :::::::::::::::::::::::::::::::::::::
            ## choose genes from pre-defined signatures
            ##

            tabPanel("Gene Signature - Kidney - Published",
              # fluidPage(
              #     column(4, actionButton("clearGeneSet", "Clear all genes")),
              #     column(4, actionButton("heatMapAction", "Plot HeatMap")),
              #  ),
               tags$hr(),
               fluidPage(
                  column(6, uiOutput("gene_signature_selection"))
               ),
               fluidPage(verbatimTextOutput('gene_set_predefined_symbols_text')),
               tags$hr(),
               fluidPage(
                  column(5, plotlyOutput("violin_gene_set_predefined"))
                  # column(7, plotlyOutput("heatmap_gene_set_predefined_centroids"))
               ),
               fluidPage(
                  # shiny::column(width = 5, plotlyOutput("violin_gene_set_predefined_clicked_gene")),
                  shiny::column(width = 12, plotlyOutput("heatmap_gene_set_predefined"))
               ),
               fluidPage(
                  shiny::column(width = 8, plotlyOutput("gene_set_predefined_rank_plot"))
               )


            ) # end tab panel





      #
      #       ## ::::::::::::::::::::::::::::::::::::::::::::
      #       ##          TAB
      #       ## ::::::::::::::::::::::::::::::::::::::::::::
      #       ##     CUSTOM GENE SIGNATURES
      #       ## ::::::::::::::::::::::::::::::::::::::::::::
      #          # Enter genes to build a custom gene signature
      #       tabPanel("Gene Signatures - Custom",
      #            # fluidPage(
      #            #     column(4, actionButton("clearGeneSet", "Clear all genes")),
      #            #     column(4, actionButton("heatMapAction", "Plot HeatMap")),
      #            #  ),
      #          tags$hr(),
      #          fluidPage(
      #             column(6, fileInput("file1", "choose text file. One gene per row.",
      #                   multiple = FALSE,
      #                   accept = c("text/csv",
      #                      "text/comma-separated-values,text/plain",
      #                      ".csv",".txt"))
      #                   )
      #                ),
      #          tags$hr(),
      #          fluidPage(
      #             verbatimTextOutput('gene_set_custom_symbols_text')
      #             ),
      #          tags$hr(),
      #          fluidPage(
      #             column(5, plotlyOutput("violin_gene_set_custom")),
      #             column(7, plotlyOutput("heatmap_gene_set_custom_centroids"))
      #          ),
      #          fluidPage(
      #             column(12, plotlyOutput("heatmap_gene_set_custom"))
      #          ),
      #          # fluidPage(
      #          #    column(5, plotlyOutput("violin_gene_set_custom_clicked_gene"))
      #          # ),
      #          tags$hr()
      #          ),
      #
      #
      #
      #
      #
      #       ## TAB 2:: GENE & CORRELATION TAB
      #       ## ---------------------------
      #       tabPanel("Gene Correlation",
      #          fluidPage(
      #             verbatimTextOutput('gene1_selection_text_t2'),
      #             verbatimTextOutput('gene2_selection_text_t2')
      #             ),
      #
      #          tags$hr(),
      #          fluidPage(
      #             column(6, plotlyOutput('violin_corr'), height=1200),
      #             column(6, plotlyOutput('scatter_gene1_gene2_corr'), height=1200)
      #             ),
      #
      #          tags$h4("Correlation table for Gene1"),
      #          fluidPage(
      #             column(6, title="Gene 1 correlated genes",DT:: dataTableOutput('corr_df'), style = "font-size:70%")
      #             ),
      #          tags$hr()
      #
      #          ),
      #
      #
      #
      #
      # ## ::::::::::::::::::::::::::::::::::::::::::::
      # ##    Pan Tissue Data
      # ## ::::::::::::::::::::::::::::::::::::::::::::
      #    tabPanel("Pan Tissue Data",
      #       fluidPage(
      #          column(12, plotlyOutput('violin_pancan'), height=1400)
      #       ),
      #       fluidPage(
      #          column(12, plotlyOutput('violin_ccle'), height=1400)
      #       )
      #    ),
      #
      # ## ::::::::::::::::::::::::::::::::::::::::::::
      # ##  TAB Transcription Factors
      # ## ::::::::::::::::::::::::::::::::::::::::::::
      #
      #       tabPanel("Transcription Factors",
      #          tags$hr(),
      #          fluidPage(
      #             column(8, plotlyOutput('violin_tf_gene'), height=1400)
      #          ),
      #          tags$hr(),
      #          fluidPage(
      #             column(4, title="TFs",DT::dataTableOutput('tf_table'), style = "font-size:70%")
      #          )
      #       ), # end lfc tab
      #
      #
      #
      # ## ::::::::::::::::::::::::::::::::::::::::::::
      # ## TAB LIMMA CONTRAST LFC TABLE
      # ## ::::::::::::::::::::::::::::::::::::::::::::
      #    tabPanel("Limma LFC table",
      #       tags$hr(),
      #       fluidPage(
      #          column(8, plotlyOutput('violin_rcc_simp_lfcgene'), height=1400)
      #       ),
      #       tags$hr(),
      #       fluidPage(
      #          column(4, title="lfc",DT::dataTableOutput('lfc_table'), style = "font-size:70%")
      #       )
      #    ), # end lfc tab
      #
      # ## ::::::::::::::::::::::::::::::::::::::::::::
      # ##    TAB LIMMA VOLCANO PLOTS
      # ## ::::::::::::::::::::::::::::::::::::::::::::
      #    tabPanel("Limma Volcanos",
      #       fluidRow(
      #            #column(4, radioButtons("pType", "Choose plot type:", choiceValues = list("A","B"),
      #            #choiceNames=list("T-statistics. Both 1st vs 2nd selected coefs", "Volcano. p-val vs FoldChange. 1st coef only"))),
      #            column(3, sliderInput(inputId = "gene_n_volcano", label = "Select n genes", min=0, max=2500, ticks = T, step = 50, value = 250)),
      #            column(3, sliderInput(inputId = "lfc_cut_volcano", label = "Log2 Fold change cut", min=0, max=6, ticks = T, step = 0.5, value = 1))
      #            ),
      #       tags$hr(),
      #
      #       # All contrasts-plots for the selected primary Coeffient
      #       tags$hr(),
      #       #verbatimTextOutput('reactive_clicked_gene_volcano_text'),
      #       fluidPage(
      #          shiny::column(width = 8, plotlyOutput("violin_rcc_simp_volcano"))
      #           ),
      #       tags$hr(),
      #       fluidPage(
      #          column(12, plotlyOutput('limma_volcano_wrap', height = 900))
      #          )
      #       ) # end tab Limma stats Plots
      #
      #
      #
      #
      #



         ) #end tabset panel
      ) # end main panel
   ) # end sidebarLayout
  ) # end UI





server <- function(input, output, session) {



## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##   SIDEBAR PANEL SELECTIONS & ASSOCIATED OBSERVE EVENTS
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   ## Selection Choices for the sidebar panel
   ## genes can also be selected & updated in some lists - 'updateVarSelectInput' associated with different 'sessions'


   ## DATA SET SELECTIION
   ## ----------------
   output$dataset_selection <- renderUI({
      selectInput(inputId = "dataset", label = "Select Dataset:",
                  choices=dataset_selection_vec, selected = dataset_selection_vec[1], multiple = F)
   })



   # ## REACTIVE - DATASET
   # reactive_eset <- reactive({
   #    eset_list[[input$dataset]]
   # })

   # reactive_gene_selection <- reactive({
   #    gene_selection_list[[input$dataset]]
   # })




   ## GENE 1 SELECTION
   ## ----------------
   output$gene1_selection <- renderUI({
            selectInput(inputId = "gene1", label = "Select Gene 1:",
                   choices=gene_selection_list[[input$dataset]], selected = gene_selection_list[[input$dataset]][1], multiple = F)
          })


   ## observe & update selected gene1 if gene sected in Limma LFC table
   observeEvent(input$lfc_table_rows_selected, {
      updateVarSelectInput(session, "gene1", selected = limmma_res_list[[input$dataset]][input$lfc_table_rows_selected,]$ENSG)
      })

   # # ## observe & update selected gene1 if gene sected in TF  table
   # # observeEvent(input$tf_table_rows_selected, {
   # #    updateVarSelectInput(session, "gene1", selected = tf_df[input$tf_table_rows_selected,]$ENSG)
   # #    })
   # #
   #
   #
   #
   # ## observe & update if clicked gene in limma volcano
   #  observeEvent(event_data("plotly_click", source = "limma_volcano_wrap"), {
   #    event_data <- event_data("plotly_click", source = "limma_volcano_wrap")
   #    updateVarSelectInput(session, "gene1", selected = event_data$key)
   # })
   # ## observe & update if clicked gene in limma limma_t_wrap
   #  observeEvent(event_data("plotly_click", source = "limma_t_wrap"), {
   #    event_data <- event_data("plotly_click", source = "limma_t_wrap")
   #    updateVarSelectInput(session, "gene1", selected = event_data$key)
   # })
   #
   #
   #
   # ## observe & update selected gene1 if gene sected in gene::mirna correlation table
   #  observeEvent(input$corr_df_mirna_rows_selected, {
   #    updateVarSelectInput(session, "gene1", selected = fData(eset_list[[input$dataset]])[input$corr_df_mirna_rows_selected,]$ENSG)
   #    })
   #
   #
   # ## observe & update gene if clicked in methylation limma table
   # ## NOTE!  If multiple genes associated to same CpG, then only change to the first of htese  genes
   #  observeEvent(input$lfc_table_me_rows_selected, {
   #       o <- input$lfc_table_me_rows_selected
   #       o_gene <- strsplit(limma_res_me_full[o,]$SYMBOL, split="[|]")[[1]][1]
   #       updateVarSelectInput(session, "gene1", selected = fData(eset_list[[input$dataset]])[which(fData(eset_list[[input$dataset]])$SYMBOL==o_gene), "ENSG"])
   #    })
   #
   #

   #
   ## GENE 2 SELECTION
   ## ----------------
   output$gene2_selection <- renderUI({
      selectInput(inputId = "gene2", label = "Select Gene 2:",
                  choices=gene_selection_list[[input$dataset]], selected = gene_selection_list[[input$dataset]][2], multiple = F)
   })
   # ## observe & update selected gene2 if gene sected in gene correlation table
   #  observeEvent(input$corr_df_rows_selected, {
   #    updateVarSelectInput(session, "gene2", selected = limmma_res_list[[input$dataset]][input$corr_df_rows_selected,]$ENSG)
   #    })
   #



   ## GENE 1 & 2 SELECTION TEXTS and HTMLs
   ## -------------------------------------
   ## gene_selection_texts
   observeEvent(input$gene1, {
      if(length(input$gene1>0)){
         output$gene1_selection_text_t1 <- output$gene1_selection_text_t2 <- output$gene1_selection_text_t3 <- renderText(paste(fData(eset_list[[input$dataset]])[input$gene1,], collapse = ", "))
      }else{output$gene1_selection_text_t1 <- NULL}
      })

   ## gene 1 HTML link
   observeEvent(input$gene1, {
      output$html_link <- renderUI({
                a("Human Protein Atlas Link", href=paste("https://www.proteinatlas.org/",input$gene1,"/tissue/kidney", sep=""), target="_blank")
        })
   })




   #
   #
   #    ## observe & update if clicked sample in scatterplot catter_annot_num_gene1
   #  observeEvent(event_data("plotly_click", source = "scatter_annot_num_gene1"), {
   #    event_data <- event_data("plotly_click", source = "scatter_annot_num_gene1")
   #    updateVarSelectInput(session, "sample_id", selected = event_data$key)
   # })
   #
   #     ## observe & update if clicked sample in scatterplot scatter_gene1_gene2
   #  observeEvent(event_data("plotly_click", source = "scatter_gene1_gene2"), {
   #    event_data <- event_data("plotly_click", source = "scatter_gene1_gene2")
   #    updateVarSelectInput(session, "sample_id", selected = event_data$key)
   # })
   #
   #     ## observe & update if clicked sample in scatterplot scatter_gene1_gene2_corr
   #  observeEvent(event_data("plotly_click", source = "scatter_gene1_gene2_corr"), {
   #    event_data <- event_data("plotly_click", source = "scatter_gene1_gene2_corr")
   #    updateVarSelectInput(session, "sample_id", selected = event_data$key)
   # })
   #
   #
   #
   #
   #
   #
   #  ## observe & update if clicked sample in scatterplot scatter_purity_methylation
   #  observeEvent(event_data("plotly_click", source = "scatter_purity_methylation"), {
   #    event_data <- event_data("plotly_click", source = "scatter_purity_methylation")
   #    updateVarSelectInput(session, "sample_id", selected = event_data$key)
   # })
   #  ## observe & update if clicked sample in scatterplot scatter_methylation_beta_vs_gex
   #  observeEvent(event_data("plotly_click", source = "scatter_methylation_beta_vs_gex"), {
   #    event_data <- event_data("plotly_click", source = "scatter_methylation_beta_vs_gex")
   #    updateVarSelectInput(session, "sample_id", selected = event_data$key)
   # })
   #
   #
   #
   #  ## observe & update if clicked sample in scatterplot scatter_gene1_gene2_selectedSample
   #  observeEvent(event_data("plotly_click", source = "scatter_gene1_gene2_selectedSample"), {
   #    event_data <- event_data("plotly_click", source = "scatter_gene1_gene2_selectedSample")
   #    updateVarSelectInput(session, "sample_id", selected = event_data$key)
   # })
   #  ## observe & update if clicked sample in scatterplot scatter_gene1_gex_cbs_selectedSample
   #  observeEvent(event_data("plotly_click", source = "scatter_gene1_gex_cbs_selectedSample"), {
   #    event_data <- event_data("plotly_click", source = "scatter_gene1_gex_cbs_selectedSample")
   #    updateVarSelectInput(session, "sample_id", selected = event_data$key)
   # })
   #
   #
   #
   #
   # ## ANNOTATION SELECTIONS (binned & numeric sample annotations)
   # ## ------------------------------------------------
   output$annotation_bin_selection <- renderUI({
            selectInput(inputId = "annot_bin", label = "Primary annotation:",
                   choices=pdata_column_types_list[[input$dataset]][['bin']], selected = "group")
          })


    # annot_num_selection
    output$annotation_num_selection <- renderUI({
            selectInput(inputId = "annot_num", label = "Primary numeric annotation:",
                        choices=pdata_column_types_list[[input$dataset]][['num']], selected = "sample_order")
          })
   #
   #
   #
   #  # ## REACTIVE - annot_bin_sub selection choices
   #  #   # Subgroup of the primary annot_bin
   #  #   annotation_bin_sub_choices <- reactive({
   #  #         as.character(unique(as.character(pdata_list[[input$dataset]][,input$annot_bin])))
   #  #      })
   #
   # # # annotation_bin_sub_selection
   # #     output$annotation_bin_sub_selection <- renderUI({
   # #          selectInput(inputId = "annot_bin_sub", label = "Primary annot - subgroup:",
   # #                 choices=annotation_bin_sub_choices(), selected = annotation_bin_sub_choices()[1])
   # #        })
   # #
   #
   #
   # ## Limma CONTRAST SELECTIONS (ALL TABS)
   # ## ------------------------------
   # # From pre-calculated Limmas - select the comparisons (contrast) to be visualized in tables and plots
   #
   #     output$contrast_1_selection <- renderUI({
   #          selectInput(inputId = "c1", label = "Limma, Coefficient 1:",
   #                 choices=contrast_1_token_list[[input$dataset]], selected = contrast_1_token_list[[input$dataset]][1])
   #        })
   #
   # output$contrast_2_selection <- renderUI({
   #          selectInput(inputId = "c2", label = "Limma, Coefficient 2:",
   #                 choices=contrast_ref_token_list[[input$dataset]][-match(input$c1, contrast_ref_token_list[[input$dataset]])], selected = contrast_ref_token_list[[input$dataset]][1])
   #        })
   #
   #
   #
   # ## REACTIVE - SELECTED CONTRASTS (ALL TABS)
   # ## ---------------------------
   # # the combination (comparison) of selected contrast 1 and contrast 2
   # reactive_contrast <- reactive({
   #    paste0(input$c1, " - " , input$c2)
   # })
   # reactive_c1_contrasts <- reactive({
   #    paste0(input$c1, " - " ,contrast_ref_token_list[[input$dataset]][!contrast_ref_token_list[[input$dataset]]==input$c1])
   # })
   #  reactive_c2_contrasts <- reactive({
   #    paste0(input$c2, " - " ,contrast_ref_token_list[[input$dataset]][!contrast_ref_token_list[[input$dataset]]==input$c2])
   # })
   # reactive_c1 <- reactive({
   #   input$c1
   # })
   #  reactive_c2 <- reactive({
   #    input$c2
   # })
   # reactive_c1c2_refs <- reactive({
   #    if(length(reactive_c2())){
   #             contrast_ref_token_list[[input$dataset]][!contrast_ref_token_list[[input$dataset]] %in% c(reactive_c1(), reactive_c2())]
   #       }else{NULL}
   #    })
   #
   # ## CONTRAST AS OUTPUT TEXT
   # ## -------------------------
   # output$contrast_selection_text1 <- output$contrast_selection_text2 <- renderText(reactive_contrast())
   #








    ## :::::::::::::::::::::::::::::::::::::
    ##      Primary Gene plot tab - RCC
    ## :::::::::::::::::::::::::::::::::::::



      ## scatter_annot_num_gene1
      ## -------------------------
         # scatter of gene 1 vs selected numeric annotation
      observe({
         if(length(input$gene1)){
            output$scatter_annot_num_gene1 <- renderPlotly({
               pdata <- pdata_list[[input$dataset]]

               p <- plotViolinScatterRankBarGeneric(
                  x = data.frame(
                     sample_id = pdata$sample_id,
                     dim1 = as.numeric(exprs(eset_list[[input$dataset]])[input$gene1, pdata$sample_id]),
                     dim2 = pdata[, input$annot_num]
                  ),
                  pdata  = pdata %>% dplyr::select(c("sample_id", !!!input$annot_bin,  !!!input$annot_num)),
                  pdata.column = input$annot_bin,
                  color.key=color_subtypes,shape.key = shape_subtypes,
                  y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL,
                  x.title = input$annot_num,
                  #plot.title = input$annot_num
                  plot.title = "Scatter: gene 1 vs Primary Numeric Annotation"
                  )

               p <- p +
                  aes(text = paste("Sample: ", sample_id, "\nGroup: ", annotation), key=sample_id)
               # ggplotly(p, tooltip = "text")
               ggplotly(p, tooltip = "text", source="scatter_annot_num_gene1") %>% layout(dragmode = "click")

            })}
      })

      ## scatter_gene1_gene2
      ## ----------------------
         # scatter of gene1 vs gene2
      observe({
         if(length(input$gene1)){
           output$scatter_gene1_gene2 <- renderPlotly({
               pdata <- pdata_list[[input$dataset]]
               print(input$gene1)
               print(input$gene2)
               str(exprs(eset_list[[input$dataset]])[input$gene1, pdata$sample_id] )
               p <- plotViolinScatterRankBarGeneric(
                  x = data.frame(
                     sample_id = pdata$sample_id,
                     dim1 = as.numeric(exprs(eset_list[[input$dataset]])[input$gene1, pdata$sample_id] ),
                     dim2 = as.numeric(exprs(eset_list[[input$dataset]])[input$gene2, pdata$sample_id] )
                     ),
                  pdata  = pdata %>% dplyr::select(c("sample_id", !!!input$annot_bin, !!!input$annot_num)),
                  pdata.column = input$annot_bin,
                  color.key=color_subtypes,shape.key = shape_subtypes,
                  y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL,
                  x.title = fData(eset_list[[input$dataset]])[input$gene2,]$SYMBOL,
                  plot.title = paste(input$annot_bin))

               p <- p +
                  aes(text = paste("Sample: ", sample_id, "\nGroup: ", annotation), key=sample_id)
               #ggplotly(p, tooltip = "text")
               ggplotly(p, tooltip = "text", source="scatter_gene1_gene2") %>% layout(dragmode = "click")

            })}
      })

      ## violin_annot_bin_gene1
      ## -------------------------
         # Violin of gene1 split into selected binned annotaton
      observe({
         if(length(input$gene1)){
            output$violin_annot_bin_gene1 <- renderPlotly({
               pdata <- pdata_list[[input$dataset]] %>% dplyr::filter(sample_id %in% pdata_list[[input$dataset]]$sample_id)
               p <- plotViolinScatterRankBarGeneric(
                  x = data.frame(
                    sample_id = pdata$sample_id,
                    dim1 = exprs(eset_list[[input$dataset]])[input$gene1,]),
                  pdata = pdata %>% dplyr::select(c("sample_id",  !!!input$annot_bin, !!!input$annot_num)),
                  pdata.column = input$annot_bin,
                  color.key=color_subtypes, y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL,
                  plot.title = "gene 1 vs Primary Annotation"
                  )

               p <- p + aes(text = paste("Group: ", annotation))
               ggplotly(p, tooltip = "text")
               })}
      }) #

   #
   #
   #
   #
   #
   #
   #
   #
   #
   #
   #
   #
   #
   #
   #
   #
   #
   # # :::::::::::::::::::::::::::::::::::::
   # #      Gene signature Explorer - Custom signaure
   # # :::::::::::::::::::::::::::::::::::::::::::::::
   # # Plot genes/scores of pre-defined & custom signatures
   #
   #
   # ## reactive gene_set_selection_file (using file)
   # ## ----------------------
   #    # Custom signature input
   #    # output$gene_set_file_selection <- renderTable({
   # reactive_gene_set_file <- reactive({
   #        # input$file1 will be NULL initially. After the user selects
   #        # and uploads a file, head of that data file by default,
   #        # or all rows if selected, will be shown.
   #       req(input$file1)
   #       df <- read.delim(input$file1$datapath, header = TRUE, sep = "\t")
   #       my_symbols <- as.character(df[,1])
   #       my_symbols <- unique(intersect(my_symbols, fData(eset_list[[input$dataset]])$SYMBOL))
   #       my_ensgs <- as.character(fData(eset_list[[input$dataset]])$ENSG[match(my_symbols, fData(eset_list[[input$dataset]])$SYMBOL)])
   #       #paste(my_ensgs)
   #       return(my_ensgs)
   #       })
   #
   #
   # ## gene_set_custom_symbols_text
   # ## -------------------------
   #    ## Generate verbatrim text output for custom signature
   # observe({
   #    my_ensgs <- reactive_gene_set_file()
   #    my_symbols <- fData(eset_list[[input$dataset]])[my_ensgs,]$SYMBOL
   #    if(length(my_symbols)){
   #       output$gene_set_custom_symbols_text <- renderText(paste0(my_symbols, collapse = ", "))
   #    }
   # })
   #
   # ## violin_gene_set_custom
   # ## --------------------
   # # Violin for custom genes set
   # observe({
   #
   #       my_ensgs <- reactive_gene_set_file()
   #
   #        if(length(my_ensgs)){
   #         output$violin_gene_set_custom <- renderPlotly({
   #             # req(input$file1)
   #             #    df <- read.delim(input$file1$datapath, header = TRUE, sep = "\t")
   #             #    # return(head(df))
   #             #    my_symbols <- as.character(df[,1])
   #             #    my_symbols <- unique(intersect(my_symbols, fData(eset_list[[input$dataset]])$SYMBOL))
   #             #    my_ensgs <- as.character(fData(eset_list[[input$dataset]])$ENSG[match(my_symbols, fData(eset_list[[input$dataset]])$SYMBOL)])
   #             #    print(my_ensgs)
   #             pdata <- pdata_list[[input$dataset]] %>% dplyr::select(c("sample_id", "tax_simp2", !!!input$annot_bin))
   #             p <- plotViolinScatterRankBarGeneric(
   #                gene.id.signature = my_ensgs,
   #                x = eset_list[[input$dataset]],
   #                pdata = pdata %>% dplyr::filter(tax_simp2 %in%
   #                       c("crtx","med","ccRCC","chRCC","pRCC_a","pRCC_b","chONC","pONC","RCC_FH","ccpRCC","tRCC","sRCC")),
   #
   #                pdata.column = input$annot_bin,
   #                color.key=color_subtypes,
   #                y.title = "signature score", plot.title = paste0("Custom gene signature (mean) vs: ", input$annot_bin)
   #                )
   #             p <- p + aes(text = paste("Group: ", annotation))
   #             ggplotly(p, tooltip = "text")
   #       })
   #       }
   #    })
   #
   # ## heatmap_gene_set_custom_centroids
   # ## --------------------------------
   # observe({
   #    my_ensgs <- reactive_gene_set_file()
   #    if(length(my_ensgs)){
   #       my_mat <- centroid_list$tax_simp3
   #       #my_mat <- cbind(my_mat, as.numeric(exprs(eset_list[[input$dataset]])[,pdata_list[[input$dataset]]$sample_id[pdata_list[[input$dataset]]$tax_simp3=="MTSCC"]]))
   #       #colnames(my_mat)[ncol(my_mat)] <- "MTSCC"
   #       head(my_mat)
   #       u <- match(levels(pdata_list[[input$dataset]]$tax_simp3), colnames(my_mat))
   #       my_mat <- my_mat[my_ensgs, u[!is.na(u)]]
   #       pdata = data.frame(
   #             sample_id = colnames(my_mat),
   #             tax_simp3 = factor(colnames(my_mat), levels=colnames(my_mat)),
   #             stringsAsFactors = F
   #              )
   #       my_mat <- t(apply(my_mat, 1, function(x) x-mean(x)))
   #       rownames(my_mat) <- featureData(eset_list[[input$dataset]][rownames(my_mat),])$SYMBOL
   #
   #    output$heatmap_gene_set_custom_centroids <- renderPlotly({
   #       p <- heatmapGeomRaster(
   #          my_mat = my_mat, pdata=pdata,
   #          plot.title = "Heatmap Custom genes. Tax Centroids",
   #          cluster.samples = F, normalize.rows = T,
   #          pdata.column.vert.lines = "tax_simp3"
   #          )
   #
   #       p <- p + aes(text = paste("Gene: ",feature_id, "\nSample: ", sample_id, "\ntax:", tax_simp3, "\n value: " , round(Z,2)), key=feature_id)
   #       ggplotly(p, tooltip = "text", source="heatmap_gene_set_custom_centroids") %>% layout(dragmode = "click")
   #    })
   #    }
   # })
   #
   #
   #
   # ## heatmap_gene_set_custom (all samples)
   # ## -----------------------
   # observe({
   #    my_ensgs <- reactive_gene_set_file()
   #    if(length(my_ensgs)){
   #     pdata = pdata_list[[input$dataset]] %>%
   #       dplyr::filter(!is.na(tax_simp2)) %>%
   #       dplyr::select(sample_id, tax, tax_simp2) %>%
   #       arrange(tax)
   #    my_es = eset_list[[input$dataset]][my_ensgs, pdata$sample_id]
   #    my_mat <- exprs(my_es) - as.numeric(apply(centroid_list_notRowCentered$tax_simp3[my_ensgs,], 1, mean))
   #    #my_mat <- t(apply(exprs(my_es), 1, function(x) x-mean(x)))
   #    rownames(my_mat) <- featureData(my_es)$SYMBOL
   #    output$heatmap_gene_set_custom <- renderPlotly({
   #       p <- heatmapGeomRaster(
   #          my_mat = my_mat, pdata=pdata,
   #          plot.title = "Heatmap Custom genes. Tax-sorted",
   #          cluster.samples = F, normalize.rows = T
   #          )
   #       p <- p + aes(text = paste("Gene: ",feature_id, "\nSample: ", sample_id, "\n value: " , round(Z,2)), key=feature_id)
   #       ggplotly(p, tooltip = "text", source="heatmap_gene_set_custom") %>% layout(dragmode = "click")
   #    })
   #    }
   # })
   #
   # # ## violin_gene_set_custom_clicked_gene
   # # ## -----------------
   # #    # Violin clicked gene , gene set published --- not vaöuable since I cannot get plotly click event to work
   # #    observe({
   # #       if(length(input$gene1)){
   # #          output$violin_gene_set_custom_clicked_gene <- renderPlotly({
   # #             pdata <- pdata_list[[input$dataset]] %>% dplyr::select(c("sample_id", "tax_simp2", !!!input$annot_bin,  !!!input$annot_num)
   # #             p <- plotViolinScatterRankBarGeneric(
   # #                gene.id = input$gene1,
   # #                x = rcc_nt,
   # #                pdata = pdata %>% dplyr::filter(tax_simp2 %in%
   # #                       c("crtx","med","ccRCC","chRCC","pRCC_a","pRCC_b","chONC","pONC","RCC_FH","ccpRCC","tRCC","sRCC")),
   # #
   # #                pdata.column = "tax_simp2",
   # #                color.key=color_subtypes,
   # #                y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL, plot.title = "GENE 1 vs Taxonomy - main groups"
   # #                )
   # #             p <- p + aes(text = paste("Group: ", annotation))
   # #             ggplotly(p, tooltip = "text")
   # #             })}
   # #    }) # end observe
   #
   #
   #


   # :::::::::::::::::::::::::::::::::::::::::::::
   #      Gene signature Explorer - Predefined
   # :::::::::::::::::::::::::::::::::::::::::::::::
      # Plot genes/scores of pre-defined signatures

      ## gene_signature_selection (genes in pre-defined signature)
      ## -------------------------
         # GENE SIGNATURE (KIDNEY PATHWAY) SELECTION
      output$gene_signature_selection <- renderUI({
               selectInput(inputId = "gene_signature", label = "Select Gene Signature:",
                      choices=names(kidney_pathways), selected = "young_litterature.Proximal", multiple = F)
             })


      ## reactive_gene_set_predefined
      ## -------------------------
         reactive_gene_set_predefined <- reactive({
            if(length(input$gene_signature)){
               # my_symbols <- as.character(kidney_pathways[[input$gene_signature]])
               # str(my_symbols)
               # my_symbols <- unique(intersect(my_symbols, fData(eset_list[[input$dataset]])$SYMBOL))
               # my_ensgs <- as.character(fData(eset_list[[input$dataset]])$ENSG[match(my_symbols, fData(eset_list[[input$dataset]])$SYMBOL)])
               my_symbols <- intersect(as.character(kidney_pathways[[input$gene_signature]]), as.character( fData(eset_list[[input$dataset]])$SYMBOL  ))
               my_ensgs <- featureNames(eset_list[[input$dataset]][fData(eset_list[[input$dataset]])$SYMBOL %in% my_symbols, ])
               str(my_ensgs)
               return(my_ensgs)
            }else{
            return()
            }
         })


      ## gene_set_predefined_symbols_text
      ## ---------------------------
      observe({
         if(length(input$gene_signature)){
            my_symbols <- intersect(as.character(kidney_pathways[[input$gene_signature]]), as.character( fData(eset_list[[input$dataset]])$SYMBOL  ))
            # my_symbols <- fData(eset_list[[input$dataset]])[my_ensgs,]$SYMBOL
            if(length(my_symbols)){
               output$gene_set_predefined_symbols_text <- renderText(paste0(my_symbols, collapse = ", "))
            }
         }
      })


      ## violin_gene_set_predefined
      ## ---------------------------------
         ## Violin for mean expression of gnes in selected pathway
      observe({
         if(length(input$gene_signature)){
            my_ensgs <- reactive_gene_set_predefined()

         if(length(my_ensgs)){
            output$violin_gene_set_predefined <- renderPlotly({
                  pdata <- pdata_list[[input$dataset]] %>% dplyr::select(c("sample_id", "tax_simp2", !!!input$annot_bin,  !!!input$annot_num))
                  p <- plotViolinScatterRankBarGeneric(
                     gene.id.signature = my_ensgs,
                     x = eset_list[[input$dataset]],
                     pdata = pdata %>% dplyr::filter(tax_simp2 %in%
                            c("crtx","med","ccRCC","chRCC","pRCC_a","pRCC_b","chONC","pONC","RCC_FH","ccpRCC","tRCC","sRCC")),

                     pdata.column = "tax_simp2",
                     color.key=color_subtypes,
                     y.title = paste0("signature score"), plot.title = paste0("SIgnature score vs main Taxonomy: ", input$gene_signature)
                     )
                  p <- p + aes(text = paste("Group: ", annotation))
                  ggplotly(p, tooltip = "text")
               })
            }}
         })

      # ## gene_set_predefined_rank_plot
      # ## ---------------------------
      #    ## gene rank plot for all genes in seleted pathway gene signature
      #  observe({
      #    if(length(reactive_c1_contrasts()) & length(input$gene_signature)){
      #       my_ensgs <- reactive_gene_set_predefined()
      #       my_symbols <- fData(eset_list[[input$dataset]])[my_ensgs,]$SYMBOL
      #          output$gene_set_predefined_rank_plot <- renderPlotly({
      #           p2 <- geneRankPlot(
      #                x = limmma_res_list[[input$dataset]] %>% dplyr::filter(contrast %in% reactive_c1_contrasts()),
      #                gene.id.column = "SYMBOL", group.column = "contrast", value.column = "t.value",
      #                n_genes = nrow(fdata_limma),
      #                color.gradient.on.rank = FALSE, ggplotply.prepare=TRUE,
      #                pathway.genes = my_symbols,
      #                pathway.name = input$gene_signature, color.key = col_key_esRankPlot)
      #
      #           p2 <- p2 + aes(text = paste("Group: ", plot.group,"\n",SYMBOL,"\nRank: ",rank.value,"\nVal: ",orig.value), key=ENSG)
      #           ggplotly(p2, tooltip = "text", source="gene_set_predefined_rank_plot") %>% layout(dragmode = "click")
      #           })
      #       }
      #  })





   # ## heatmap_gene_set_predefined_centroids
   # ## --------------------------------
   # observe({
   #    if(length(input$gene_signature)){
   #       my_ensgs <- reactive_gene_set_predefined()
   #       if(length(my_ensgs)){
   #          my_mat <- centroid_list$tax_simp3
   #          #my_mat <- cbind(my_mat, as.numeric(exprs(rcc_nt[,pdata_list[[input$dataset]]$sample_id[pdata_list[[input$dataset]]$tax_simp3=="MTSCC"]])))
   #          #colnames(my_mat)[ncol(my_mat)] <- "MTSCC"
   #          head(my_mat)
   #          u <- match(levels(pdata_list[[input$dataset]]$tax_simp3), colnames(my_mat))
   #          my_mat <- my_mat[my_ensgs, u[!is.na(u)]]
   #          pdata = data.frame(
   #                sample_id = colnames(my_mat),
   #                tax_simp3 = factor(colnames(my_mat), levels=colnames(my_mat)),
   #                stringsAsFactors = F
   #                 )
   #          my_mat <- t(apply(my_mat, 1, function(x) x-mean(x)))
   #          rownames(my_mat) <- featureData(eset_list[[input$dataset]][rownames(my_mat),])$SYMBOL
   #
   #       output$heatmap_gene_set_predefined_centroids <- renderPlotly({
   #          p <- heatmapGeomRaster(
   #             my_mat = my_mat, pdata=pdata,
   #             plot.title = "Heatmap Custom genes. Taxonomy centroids",
   #             cluster.samples = F, normalize.rows = T,
   #             pdata.column.vert.lines = "tax_simp3"
   #             )
   #          p <- p + aes(text = paste("Gene: ",feature_id, "\nSample: ", sample_id, "\ntax:",tax_simp3, "\n value: " , round(Z,2)), key=feature_id)
   #          ggplotly(p, tooltip = "text", source="heatmap_gene_set_predefined_centroids") %>% layout(dragmode = "click")
   #       })
   #       }}
   # })
   #
   #
   #


   # ## heatmap_gene_set_predefined
   # ## -------------------------------
   # observe({
   #    if(length(input$gene_signature)){
   #       my_ensgs <- reactive_gene_set_predefined()
   #       if(length(my_ensgs)){
   #          pdata = pdata_list[[input$dataset]] %>%
   #             dplyr::filter(!is.na(tax_simp2)) %>%
   #             dplyr::select(sample_id, tax, tax_simp2) %>%
   #             arrange(tax)
   #          my_es = eset_list[[input$dataset]][my_ensgs, pdata$sample_id]
   #          my_mat <- exprs(my_es) - as.numeric(apply(centroid_list_notRowCentered$tax_simp3[my_ensgs,], 1, mean))
   #          #my_mat <- t(apply(exprs(my_es), 1, function(x) x-mean(x)))
   #          rownames(my_mat) <- featureData(my_es)$SYMBOL
   #
   #          output$heatmap_gene_set_predefined <- renderPlotly({
   #             p <- heatmapGeomRaster(
   #                my_mat = my_mat, pdata=pdata,
   #                plot.title = "Heatmap Custom genes. Tax-sorted",
   #                cluster.samples = F, normalize.rows = T
   #                )
   #             p <- p + aes(text = paste("Gene: ",feature_id, "\nSample: ", sample_id, "\n value: " , round(Z,2)), key=feature_id)
   #             ggplotly(p, tooltip = "text", source="heatmap_gene_set_predefined") %>% layout(dragmode = "click")
   #       })
   #    }}
   # })
   #




   #
   #
   #
   #
   #
   #
   # ## :::::::::::::::::::::::::::::::::::::
   # ##      Gene Correlation analysis & plots
   # ## :::::::::::::::::::::::::::::::::::::
   #    # For correlation analysis - sample groups (taonomies) can be added/removed to perform correlation on a subset
   #
   #    ## es_reactive
   #    ## ----------------------------------
   #
   #
   #    ## corr_df
   #    ## download_corr_df
   #    ## -------------------
   #       # Correlation table for selected gene vs within taxonomies. A download hanlde table for download of correlations
   #    observe({
   #       if(length(input$gene2)>0){
   #          cor_tab_rcc <- fData(eset_list[[input$dataset]])
   #          cor_tab_rcc$r <- suppressWarnings(round(as.numeric(cor(as.numeric(exprs(es_reactive()[input$gene1,])), t(exprs(es_reactive())))), 3))
   #
   #          output$corr_df <- DT::renderDataTable(cor_tab_rcc, server = TRUE, selection = 'single', options = list(pageLength = 25))
   #
   #          output$download_corr_df <- downloadHandler(
   #              filename = function() {
   #                paste("cor_table_", gene_table$SYMBOL[input$gene1], ".txt", sep = "", quote=F)
   #              },
   #              content = function(file) {
   #                write.table(cor_tab_rcc[order(cor_tab_rcc$r, decreasing = T),], file, sep="\t", row.names = FALSE)
   #             })
   #       }
   #    })
   #
   #    ## violin_corr
   #    ## -----------------
   #       # viloin of gene selected from the correlation table
   #    observe({
   #       if(length(input$gene2)){
   #          output$violin_corr <- renderPlotly({
   #             pdata <- pdata_list[[input$dataset]]
   #             p <- plotViolinScatterRankBarGeneric(
   #                gene.id = input$gene2,
   #                x = eset_list[[input$dataset]],
   #                pdata = pdata %>% dplyr::filter(tax_simp2 %in%
   #                       c("crtx","med","ccRCC","chRCC","pRCC_a","pRCC_b","chONC","pONC","RCC_FH","ccpRCC","tRCC","sRCC")) %>%
   #                         dplyr::select(c("sample_id", "tax_simp2", !!!input$annot_bin,  !!!input$annot_num)),
   #                pdata.column = "tax_simp2",
   #                color.key=color_subtypes,
   #                y.title = fData(eset_list[[input$dataset]])[input$gene2,]$SYMBOL,
   #                plot.title = paste("RCC Taxonomy Simplified")
   #                )
   #             p <- p + aes(text = paste("Group: ", annotation))
   #             ggplotly(p, tooltip = "text")
   #             })}
   #    }) # end observe
   #
   #    ## scatter_gene1_gene2_corr
   #    ## ------------------------
   #       # Plot scatter gene1 vs gene correlated
   #    observe({
   #       if(length(input$gene2)){
   #
   #          output$scatter_gene1_gene2_corr <- renderPlotly({
   #             pdata <- pdata_list[[input$dataset]]
   #             p <- plotViolinScatterRankBarGeneric(
   #                x = data.frame(
   #                   sample_id = sampleNames(es_reactive()),
   #                   dim1 = as.numeric(exprs(es_reactive()[input$gene1, ])),
   #                   dim2 = as.numeric(exprs(es_reactive()[input$gene2, ]))
   #                   ),
   #                pdata = pdata, pdata.column = input$annot_bin,
   #                color.key=color_subtypes,
   #                shape.key = shape_subtypes,
   #                x.title = fData(eset_list[[input$dataset]])[input$gene2,]$SYMBOL,
   #                y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL,
   #                plot.title = "Gene 1 vs Correlated gene")
   #
   #             p <- p +
   #                aes(text = paste("Sample: ", sample_id, "\ntax_simp:",tax_simp2, "\nGroup: ", annotation), key=sample_id)
   #             ggplotly(p, tooltip = "text", source="scatter_gene1_gene2_corr") %>% layout(dragmode = "click")
   #             #ggplotly(p)
   #          })}
   #    })
   #
   #
   #
   #
   #
   #
   #
   # ## :::::::::::::::::::::::::::::::::::::
   # ##      Limma volcano plots
   # ## :::::::::::::::::::::::::::::::::::::
   #    ## LIMMA VOLCANO PLOTS (Volcano)
   #
   #    ## limma_volcano_wrap
   #    ## ----------------------
   #       ## wrapped volcano plots (facet_wrap) of pre-calculated Limma results- Log2 diff vs p-value
   #    observe({
   #       if(length(reactive_c1_contrasts())){
   #          output$limma_volcano_wrap <- renderPlotly({
   #
   #             p4 <- plotViolinScatterRankBarGeneric(
   #               x = data.frame(limmma_res_list[[input$dataset]] %>%
   #                     dplyr::filter(contrast %in% reactive_c1_contrasts()) %>%
   #                     dplyr::filter(abs(lfc)>= input$lfc_cut_volcano) %>%
   #                     dplyr::filter(gene_padj_BH_call !="ns") %>%
   #                     dplyr::mutate(dim1 = -1*log10(p.value), dim2=lfc, annotation = gene_padj_BH_call) %>%
   #                     group_by(contrast) %>%
   #                     top_n(n = input$gene_n_volcano, wt = -p.value)),
   #                ylim.rescale=c(0,NA),
   #                color.key = c(Downreg="deepskyblue2", Upreg="indianred"),
   #                x.title = "Log2 fold diff. expression", y.title = "-log10 p.value",
   #                zero.lines=T
   #             )
   #             p4 <- p4 +
   #                facet_wrap(~contrast) + aes(key=ENSG) +
   #                aes(text = paste("Contrast: ", contrast,"\n", ENSG, "\n", SYMBOL, "\np.val: ", signif(p.value, 2)))
   #             ggplotly(p4, tooltip = "text", source="limma_volcano_wrap") %>% layout(dragmode = "click")
   #          })
   #       }
   #    })
   #
   #
   #    ## violin_rcc_simp_volcano
   #    ## ---------------------------
   #       # Observe reactive click event for gene in rankplot - then create violin
   #    observe({
   #          if(length(input$gene1)){
   #             output$violin_rcc_simp_volcano <- renderPlotly({
   #                pdata <- pdata_list[[input$dataset]] %>% dplyr::select(c("sample_id", "tax_simp2", !!!input$annot_bin,  !!!input$annot_num))
   #                p <- plotViolinScatterRankBarGeneric(
   #                   gene.id = input$gene1,
   #                   x = eset_list[[input$dataset]],
   #                   pdata = pdata %>% dplyr::filter(tax_simp2 %in%
   #                          c("crtx","med","ccRCC","chRCC","pRCC_a","pRCC_b","chONC","pONC","RCC_FH","ccpRCC","tRCC","sRCC")),
   #
   #                   pdata.column = "tax_simp2",
   #                   color.key=color_subtypes,
   #                   y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL, plot.title = "GENE 1 vs Taxonomy - main groups"
   #                   )
   #                p <- p + aes(text = paste("Group: ", annotation))
   #                ggplotly(p, tooltip = "text")
   #                })}
   #       }) # end observe
   #
   #
   # ## :::::::::::::::::::::::::::::::::::::
   # ##      Limma analyses - T-scatter plots
   # ## :::::::::::::::::::::::::::::::::::::
   #    # LIMMA T-SCATTER PLOTS
   #
   #    ## limma_t_wrap
   #    ## ----------------
   #    observe({
   #       if(!is.null(reactive_c1c2_refs) & !any(reactive_c2() %in% c("crtx","med")) ){
   #          #c1 <- "ccpRCC"
   #          #c2 <- "chRCC"
   #
   #             print(reactive_c1c2_refs())
   #             my_c1s <- paste0(reactive_c1(), " - ", reactive_c1c2_refs())
   #             my_c2s <- paste0(reactive_c2(), " - ", reactive_c1c2_refs())
   #             #print(my_c2s)
   #
   #             y1 <- data.frame(limmma_res_list[[input$dataset]] %>%
   #                dplyr::filter(contrast %in% my_c1s) %>%
   #                mutate(ENSG = as.character(ENSG)) %>%
   #                mutate(sample_id=paste0(c2,ENSG)) %>%
   #                mutate(dim1 = t.value, annotation = gene_padj_BH_call) %>%
   #                dplyr::select(sample_id, dim1, c1, c2, contrast, lfc, ENSG, SYMBOL, annotation) %>%
   #                group_by(contrast) %>% top_n(n = input$gene_n_t, wt = abs(dim1))  ) %>% droplevels()
   #             #str(y1)
   #             #print(head(y1))
   #
   #             y2 <- limmma_res_list[[input$dataset]] %>%
   #                dplyr::filter(contrast %in% as.character(my_c2s)) %>%
   #                mutate(ENSG = as.character(ENSG)) %>%
   #                dplyr::filter(ENSG %in% as.character(y1$ENSG)) %>%
   #                mutate(dim2 = t.value) %>%
   #                mutate(sample_id=paste0(c2, ENSG)) %>% droplevels()
   #
   #             #print(head(y2))
   #
   #             y <- dplyr::left_join(y1, y2 %>% dplyr::select(sample_id, dim2, c1, contrast), by="sample_id")
   #             #str(y)
   #             rm(y1,y2)
   #
   #             my_ylim.rescale <- c(-50,50)
   #             my_xlim.rescale <- c(-50,50)
   #
   #             if(nrow(y)>=1){
   #             output$limma_t_wrap <- renderPlotly({
   #                print(head(y))
   #                p <- plotViolinScatterRankBarGeneric(
   #                   x =  data.frame(y),
   #                   ylim.rescale=my_ylim.rescale, xlim.rescale=my_xlim.rescale,
   #                   color.key = c(Downreg="deepskyblue2", Upreg="indianred"),
   #                   # x.title = "", y.title = "",
   #                   zero.lines=T,
   #                   x.title = input$c2, y.title = input$c1
   #                   )
   #                p <- p +
   #                   facet_wrap(~c2) + aes(key=ENSG) +
   #                   aes(text = paste("t.values for:\n",contrast.x, " vs ", contrast.y, "\n", ENSG, "\n", SYMBOL))
   #                ggplotly(p, tooltip = "text", source="limma_t_wrap") %>% layout(dragmode = "click")
   #             })
   #             }
   #          }
   #       })
   #
   #
   #    ## violin_rcc_simp_t
   #    ## ----------------------------------------------------
   #       # Observe reactive click event for rankplot - then create violin
   #    observe({
   #          if(length(input$gene1)){
   #             output$violin_rcc_simp_t <- renderPlotly({
   #                pdata <- pdata_list[[input$dataset]] %>% dplyr::select(c("sample_id", "tax_simp2", !!!input$annot_bin,  !!!input$annot_num))
   #                p <- plotViolinScatterRankBarGeneric(
   #                   gene.id = input$gene1,
   #                   x = eset_list[[input$dataset]],
   #                   pdata = pdata %>% dplyr::filter(tax_simp2 %in%
   #                          c("crtx","med","ccRCC","chRCC","pRCC_a","pRCC_b","chONC","pONC","RCC_FH","ccpRCC","tRCC","sRCC")),
   #
   #                   pdata.column = "tax_simp2",
   #                   color.key=color_subtypes,
   #                   y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL, plot.title = "GENE 1 vs Taxonomy - main groups"
   #                   )
   #                p <- p + aes(text = paste("Group: ", annotation))
   #                ggplotly(p, tooltip = "text")
   #                })}
   #       }) # end observe
   #
   #
   #
   #
   # ## :::::::::::::::::::::::::::::::::::::
   # ##      Limma analyses - Log fold difference tables
   # ## :::::::::::::::::::::::::::::::::::::
   #    # between selected limma contrasts in sidebar menu
   #
   #  ## lfc_table
   #  ## --------------------
   #     ## LFC table - c1 vs c2
   #     observe({
   #       if(length(reactive_c2())){
   #          print(reactive_contrast())
   #          output$lfc_table <- DT::renderDataTable(
   #               data.frame(limmma_res_list[[input$dataset]] %>%
   #                     dplyr::filter(as.character(contrast)==reactive_contrast()) %>% droplevels() %>%
   #                     dplyr::select(contrast, ENSG, SYMBOL, lfc, gene_padj_BH_call, p.value, t.value) %>%
   #                     mutate(lfc = signif(lfc, 3), p.value = signif(p.value, 3), t.value = signif(t.value, 3))
   #                  ), options = list(pageLength = 25), server = TRUE, selection = 'single')
   #       }
   #       # output$download_lfc_table <- downloadHandler(
   #       #           filename = function() {
   #       #             paste("lfc_table_", c1,"_vs_",c2, ".txt", sep = "")
   #       #           },
   #       #           content = function(file) {
   #       #             write.table(data.frame(cbind(t_df[,c("ENSG","SYMBOL")], round(centroid_list$tax_simp2[,c(c1,c2)], 2)) %>% mutate(lfc_diff = .[[3]] - .[[4]])), file, quote=F, sep="\t", row.names = FALSE)
   #       #          })
   #     })
   #
   #    ## violin_rcc_simp_lfcgene
   #    ## -----------------------
   #     observe({
   #          if(length(input$gene1)){
   #             output$violin_rcc_simp_lfcgene <- renderPlotly({
   #                pdata <- pdata_list[[input$dataset]] %>% dplyr::select(c("sample_id", "tax_simp2", !!!input$annot_bin,  !!!input$annot_num))
   #                p <- plotViolinScatterRankBarGeneric(
   #                   gene.id = input$gene1,
   #                   x = eset_list[[input$dataset]],
   #                   pdata = pdata %>% dplyr::filter(tax_simp2 %in%
   #                          c("crtx","med","ccRCC","chRCC","pRCC_a","pRCC_b","chONC","pONC","RCC_FH","ccpRCC","tRCC","sRCC")),
   #
   #                   pdata.column = "tax_simp2",
   #                   color.key=color_subtypes,
   #                   y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL, plot.title = "GENE 1 vs Taxonomy - main groups"
   #                   )
   #                p <- p + aes(text = paste("Group: ", annotation))
   #                ggplotly(p, tooltip = "text")
   #                })}
   #       }) # end observe
   #
   #
   #
   # ## :::::::::::::::::::::::::::::::::::::
   # ##      Transcription factors
   # ## :::::::::::::::::::::::::::::::::::::
   #
   #    ## tf_table
   #    ## --------------------
   #     output$tf_table <- DT::renderDataTable(
   #               data.frame(tf_df), options = list(pageLength = 25), server = TRUE, selection = 'single'
   #        )
   #
   #    ## violin_tf_gene
   #    ## -----------------------
   #      observe({
   #          if(length(input$gene1)){
   #             output$violin_tf_gene <- renderPlotly({
   #                pdata <- pdata_list[[input$dataset]] %>% dplyr::select(c("sample_id", "tax_simp2", !!!input$annot_bin,  !!!input$annot_num))
   #                p <- plotViolinScatterRankBarGeneric(
   #                   gene.id = input$gene1,
   #                   x = eset_list[[input$dataset]],
   #                   pdata = pdata %>% dplyr::filter(tax_simp2 %in%
   #                          c("crtx","med","ccRCC","chRCC","pRCC_a","pRCC_b","chONC","pONC","RCC_FH","ccpRCC","tRCC","sRCC")),
   #
   #                   pdata.column = "tax_simp2",
   #                   color.key=color_subtypes,
   #                   y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL, plot.title = "GENE 1 vs Taxonomy - main groups"
   #                   )
   #                p <- p + aes(text = paste("Group: ", annotation))
   #                ggplotly(p, tooltip = "text")
   #                })}
   #       }) # end observe
   #
   #
   #
   # ## :::::::::::::::::::::::::::::::::::::
   # ##      Pan Cancer & Tissue data
   # ## :::::::::::::::::::::::::::::::::::::
   #
   #  ## violin_pancan
   #  ## --------------
   #   observe({
   #       if(length(input$gene1)){
   #          output$violin_pancan <- renderPlotly({
   #             pdata <- pdata_panCanGex %>% dplyr::select(c("sample_id", "sample_type2", !!!input$annot_bin,  !!!input$annot_num))
   #             p <- plotViolinScatterRankBarGeneric(
   #                gene.id = input$gene1,
   #                x = pan_gdc,
   #                pdata = pdata,
   #                pdata.column = "sample_type2",
   #                color.key=color_subtypes,
   #                y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL, plot.title = "GENE 1 vs TCGA PanCan sample types"
   #                )
   #             p <- p + aes(text = paste("Group: ", annotation))
   #             ggplotly(p, tooltip = "text")
   #             })}
   #    }) # end observe
   #
   #  ## violin_ccle
   #  ## --------------
   #   observe({
   #       if(length(input$gene1)){
   #          output$violin_ccle <- renderPlotly({
   #             pdata <- pdata_ccle
   #             p <- plotViolinScatterRankBarGeneric(
   #                gene.id = input$gene1,
   #                x = ccle,
   #                pdata = pdata,
   #                pdata.column = "tissue",
   #                color.key=color_subtypes,
   #                y.title = fData(eset_list[[input$dataset]])[input$gene1,]$SYMBOL, plot.title = "GENE 1 vs CCLE cell types"
   #                )
   #             p <- p + aes(text = paste("Group: ", annotation))
   #             ggplotly(p, tooltip = "text")
   #             })}
   #    }) # end observe
   #





   } # end server


# Run the application
shinyApp(ui = ui, server = server)
