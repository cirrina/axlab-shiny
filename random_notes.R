
## ADD datasets
dataset_selection_vec


///



violin_rcc_simp_gene1 // violin_gene1
violin_rcc_simp_gene2  // violin_gene2


violin_rcc_corr -- violin_corr

keep
scatter_annot_num_gene1
scatter_gene1_gene2
violin_annot_bin_gene1

barplot_annot_bin


gene_signature_selection
gene_set_predefined_symbols_text
violin_gene_set_predefined
heatmap_gene_set_predefined_centroids
heatmap_gene_set_predefined
gene_set_predefined_rank_plot


scatter_gene1_gene2_corr


REMOVE
violin_annot_bin_num
scatter_annot_num_num2
violin_annot_bin_bin_wrap


taxDataToken


mutation_tab_gene1

—Copy Number Plots - Gene
violin_rcc_simp_gene1_cbs
barplot_rcc_annot_gene1_ascn_cn
barplot_rcc_annot_gene1_ascn_cn_minor
stripchart_gene1_gex_ascn_all
stripchart_gene1_gex_ascn_minor
violin_rcc_annot_bin_gene1_cbs


sample_selection_text_t1
gene1_selection_text_t3
tsne_gex_singleSample
tsne_me_singleSample
tsne_mir_singleSample
tsne_rppa_singleSample
umap_gex_singleSample
umap_me_singleSample
umap_mir_singleSample
umap_rppa_singleSample
genome_plot
genome_plot_circle
scatter_gene1_gene2_selectedSample
scatter_gene1_gex_cbs_selectedSample
scatter_gene1_cbs_purity_selectedSample
mutation_tab_cosmic
mutation_tab


Annotation subgroup
annotation_bin_sub_selection
mutation_tab_annot_sub_summary
mutation_tab_annot_sub


TAB: FGSEA Signatures & Enrichment Tables
ea_term_plot **
   heatmap_ea
.. etc ea

TSNE
tsne_gex

….

violin_rcc_simp_t
limma_t_wrap


Methylation limma
violin_methylation, scatter_methylation_beta_vs_gex … tex


mirNA
violin_pancan_mir, violin_miRNA, scatter_gene1_mirna_corr


RPPA limma
miRNA
