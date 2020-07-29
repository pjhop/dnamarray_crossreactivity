## Format and save figures for manuscript

## MAIN FIGURES
## Figure 1 - qqplot + manhattan
qqplot <- readRDS("../data/output/figs/qqplot_mlma_loco.rds")
manhattan <- readRDS("../data/output/figs/manhattan_mlma_loco.rds")

figure1 <- cowplot::plot_grid(manhattan, qqplot, nrow = 1,
                              labels = c("A", "B"), align = "vh",
                              label_size = 10)
ggplot2::ggsave(figure1, filename = "../data/output/figs/ms/figure1.png",
                dpi = 600, width = 2*3.0, height = 3.0 * 3/4)

## Figure 3 - C9 match lengths
plot_typeI_nomismatch_noindel_correlated_methylated <- readRDS("../data/output/figs/c9_overlaps/typeI_nomismatch_noindel_correlated_methylated.rds")
plot_typeI_nomismatch_noindel_correlated_methylated_oobpvalues <- readRDS("../data/output/figs/c9_overlaps/typeI_nomismatch_noindel_correlated_methylated_oobpvalues.rds")
diff_plot_nomismatch_mismatch_indel <- readRDS("../data/output/figs/c9_overlaps/diff_plot_nomismatch_mismatch_indel.rds")
plot_typeI_mismatch_indel_correlated_methylated_oobpvalues <- readRDS("../data/output/figs/c9_overlaps/typeI_mismatch_indel_correlated_methylated_oobpvalues.rds")

figure3 <- cowplot::plot_grid(plot_typeI_nomismatch_noindel_correlated_methylated, plot_typeI_nomismatch_noindel_correlated_methylated_oobpvalues,
                              diff_plot_nomismatch_mismatch_indel,
                              plot_typeI_mismatch_indel_correlated_methylated_oobpvalues,
                              labels = c("A", "B", "C", "D"), nrow = 2, label_size = 10, align = "vh")
ggplot2::ggsave(figure3, dpi = 600, width = 3.0 * 2, height = 3.0 * 3/4 * 2,
                filename = "../data/output/figs/ms/figure3.png")


## Figure 4 - Intensity EWASs
compare_mlma_total <- readRDS("../data/output/figs/compare_EWAS_beta_EWAS_total.rds")
compare_mlma_total_oob <- readRDS("../data/output/figs/compare_EWAS_beta_EWAS_total_oob.rds")

figure4 <- cowplot::plot_grid(compare_mlma_total, compare_mlma_total_oob,
                              labels = c("A", "B"), nrow = 1, label_size = 10, align = "vh")
ggplot2::ggsave(figure4, dpi = 600, width = 3.0 * 2, height = 3.0,
                filename = "../data/output/figs/ms/figure4.png")

## Figure 5 - EPIC replication
compare_450k_epic <- readRDS("../data/output/figs/compare_450k_epic.rds")
compare_450k_epic_oob <- readRDS("../data/output/figs/compare_450k_epic_oob.rds")

figure5 <- cowplot::plot_grid(compare_450k_epic, compare_450k_epic_oob,
                              labels = c("A", "B"), nrow = 1, label_size = 10, align = "vh")
ggplot2::ggsave(figure5, dpi = 600, width = 3.0 * 2, height = 3.0,
                filename = "../data/output/figs/ms/figure5.png")


## SUPPLEMENTAL FIGURES

# Sensitivity analyses EWAS
compare_rmPCAoutliers <- readRDS("../data/output/figs/compare_rmPCAoutliers.rds")
compare_linear <- readRDS("../data/output/figs/compare_linear.rds")
compare_meta <- readRDS("../data/output/figs/compare_meta.rds")
compare_mvalues <- readRDS("../data/output/figs/compare_mvalues.rds")

sensitivity_analyses <- cowplot::plot_grid(compare_rmPCAoutliers, compare_meta,
                                   compare_linear, compare_mvalues,
                              labels = c("A", "B", "C", "D"), nrow = 2, label_size = 10, align = "vh")

ggplot2::ggsave(sensitivity_analyses, dpi = 600, width = 3.0 * 2, height = 3.0 * 2,
                filename = "../data/output/figs/ms/sensitivity_analyses.png")

# Manhattan + qqplot OOB
qqplot_oob <- readRDS("../data/output/figs/qqplot_oob.rds")
manhattan_oob <- readRDS("../data/output/figs/manhattan_oob.rds")

manhattan_qqplot_oob <- cowplot::plot_grid(manhattan_oob, qqplot_oob, nrow = 1,
                              labels = c("A", "B"), align = "vh",
                              label_size = 10)
ggplot2::ggsave(manhattan_qqplot_oob, filename = "../data/output/figs/ms/manhattan_qqplot_oob.png",
                dpi = 600, width = 2*3.0, height = 3.0 * 3/4)

# T-statistics
plot_typeI_mismatch_indel_methylated_tstat <- readRDS("../data/output/figs/c9_overlaps/typeI_mismatch_indel_methylated_tstat.rds")
ggsave(plot_typeI_mismatch_indel_methylated_tstat,
      dpi = 600,width = 3.0, height = 3.0 * 3/4, filename = "../data/output/figs/ms/tstats_typeI.png")

#  Type II plots
plot_typeII_nomismatch_noindel_correlated_methylated <- readRDS("../data/output/figs/c9_overlaps/typeII_nomismatch_noindel_correlated_methylated.rds")
plot_typeII_mismatch_indel_correlated_methylated <- readRDS("../data/output/figs/c9_overlaps/typeII_mismatch_indel_correlated_methylated.rds")

typeII_plots <- cowplot::plot_grid(plot_typeII_nomismatch_noindel_correlated_methylated,
                                   plot_typeII_mismatch_indel_correlated_methylated,
                                   nrow = 1,
                                   labels = c("A", "B"), align = "vh",
                                   label_size = 10)
ggplot2::ggsave(typeII_plots, filename = "../data/output/figs/ms/plots_c9_overlaps_typeII.png",
                dpi = 600, width = 2*3.0, height = 3.0 * 3/4)

# Manhattan plots, highlighting cross-reactive probes
manhattan_hlightcrossreactive <- readRDS("../data/output/figs/manhattan_hlightcrossreactive.rds")
manhattan_oob_hlightcrossreactive <- readRDS("../data/output/figs/manhattan_oob_hlightcrossreactive.rds")

manhattan_hlight_ib_oob <- cowplot::plot_grid(manhattan_hlightcrossreactive,
                                   manhattan_oob_hlightcrossreactive,
                                    nrow = 1,
                                    labels = c("A", "B"), align = "vh",
                                    label_size = 10)
ggplot2::ggsave(manhattan_hlight_ib_oob, filename = "../data/output/figs/ms/manhattan_hlight_ib_oob.png",
                dpi = 600, width = 2*3.0, height = 3.0 * 3/4)

# OR plot
OR_plot <- readRDS("../data/output/figs/ORs_plot.rds")

ggplot2::ggsave(OR_plot, filename = "../data/output/figs/ms/OR_plot.png",
                dpi = 600, width = 2*3.0, height = 3.0 * 3/4)


### Sensitivity ploys

# Methylated vs unmethylated type I + type II

diff_plot_unmethylated_methylated_typeI <- readRDS("../data/output/figs/c9_overlaps/diff_plot_unmethylated_methylated_typeI.rds")
diff_plot_unmethylated_methylated_typeII <- readRDS("../data/output/figs/c9_overlaps/diff_plot_unmethylated_methylated_typeII.rds")

diff_plots_U_M <- cowplot::plot_grid(diff_plot_unmethylated_methylated_typeI,
                                   diff_plot_unmethylated_methylated_typeII,
                                   nrow = 1,
                                   labels = NULL, align = "vh",
                                   label_size = 10,
                                   hjust = 0.5)
ggplot2::ggsave(diff_plots_U_M, filename = "../data/output/figs/ms/diff_plots_U_M.png",
                dpi = 600, width = 2*3.0, height = 3.0 * 3/4)

# YpG
# Type I
diff_plot_methylated_YpG_typeI <- readRDS("../data/output/figs/c9_overlaps/diff_plot_methylated_YpG_mismatch_indel_typeI.rds")
diff_plot_unmethylated_YpG_typeI <- readRDS("../data/output/figs/c9_overlaps/diff_plot_unmethylated_YpG_mismatch_indel_typeI.rds")

diff_plots_YpG_typeI <- cowplot::plot_grid(diff_plot_methylated_YpG_typeI,
                                     diff_plot_unmethylated_YpG_typeI,
                                     nrow = 1,
                                     labels = NULL, align = "vh",
                                     label_size = 10)
ggplot2::ggsave(diff_plots_YpG_typeI, filename = "../data/output/figs/ms/diff_plots_YpG_typeI.png",
                dpi = 600, width = 2*3.0, height = 3.0 * 3/4)


# Type II
diff_plot_methylated_YpG_typeII <- readRDS("../data/output/figs/c9_overlaps/diff_plot_methylated_YpG_mismatch_indel_typeII.rds")
diff_plot_unmethylated_YpG_typeII <- readRDS("../data/output/figs/c9_overlaps/diff_plot_unmethylated_YpG_mismatch_indel_typeII.rds")

diff_plots_YpG_typeII <- cowplot::plot_grid(diff_plot_methylated_YpG_typeII,
                                            diff_plot_unmethylated_YpG_typeII,
                                           nrow = 1,
                                           labels = NULL, align = "vh",
                                           label_size = 10)
ggplot2::ggsave(diff_plots_YpG_typeII, filename = "../data/output/figs/ms/diff_plots_YpG_typeII.png",
                dpi = 600, width = 2*3.0, height = 3.0 * 3/4)

## Include flanking regions
diff_plot_flanking_all <- readRDS("../data/output/figs/c9_overlaps/diff_plot_all_mismatch_indel_includeflanking.rds")
ggplot2::ggsave(diff_plot_flanking_all, filename = "../data/output/figs/ms/diff_plot_flanking_all.png",
                dpi = 600, width = 3.0, height = 3.0 * 3/4)


## mismatch/INDEL >3bp
diff_plot_3bp_all <- readRDS("../data/output/figs/c9_overlaps/diff_plot_all_mismatch_indel_3bp.rds")
ggplot2::ggsave(diff_plot_3bp_all, filename = "../data/output/figs/ms/diff_plot_3bp_all.png",
                dpi = 600, width = 3.0, height = 3.0 * 3/4)

## EPIC
manhattan_EPIC <- readRDS("../data/output/figs/manhattan_EPIC.rds")
qqplot_EPIC <- readRDS("../data/output/figs/qqplot_EPIC.rds")

manhattan_qqplot_EPIC <- cowplot::plot_grid(manhattan_EPIC, qqplot_EPIC, nrow = 1,
                              labels = c("A", "B"), align = "vh",
                              label_size = 10)

ggplot2::ggsave(manhattan_qqplot_EPIC, filename = "../data/output/figs/ms/manhattan_qqplot_EPIC.png",
                dpi = 600, width = 2*3.0, height = 3.0 * 3/4)

## Beta-plots
beta_plot_ib <- readRDS("../data/output/figs/beta_plot_ib.rds")
beta_plot_oob <- readRDS("../data/output/figs/beta_plot_oob.rds")


beta_plots <- cowplot::plot_grid(beta_plot_ib, beta_plot_oob,
                              labels = c("A", "B"), nrow = 1, label_size = 10, align = "vh")
ggplot2::ggsave(beta_plots, dpi = 600, width = 4 * 2, height = 3.0,
                filename = "../data/output/figs/ms/beta_plots.png")
