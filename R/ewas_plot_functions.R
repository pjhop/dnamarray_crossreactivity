# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## Adapted from: https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function#script
gg.manhattan <- function(stats, array = "450k", annotate_stats = FALSE, threshold = NULL, sugg = NULL,
                         hlight = NULL, annotate = FALSE,
                         title = NULL, col = NULL, ylims = NULL, text_size = 9, point_size = 1){

  if(is.null(threshold)) {
    threshold <- 0.05 / nrow(stats) # Use bonferroni by default
  }
  if(is.null(title)) {
    title <- ""
  }
  if(is.null(ylims)) {
    logp <- -log10(stats$p)
    max_logp <- max(logp)
    ylims <- c(0, max_logp + 0.15 * max_logp)
  }
  if(is.null(col)) {
    col <- c("#4EAFAF", "#2C9696", "#0F8F8F", "#057272", "#005A5A") # m
  }


  # Add annotation, add EPIC later!!
  if(!is.null(annotate_stats)) {
    ## Add annotation
    stats <- stats[,!colnames(stats) %in% c("chr", "pos", "Name", "UCSC_RefGene_Name", "bp", "gene")]
    stats <- stats %>% dplyr::left_join(annotate_stats[,c("chr", "pos", "Name", "UCSC_RefGene_Name")], by = c("Probe" = "Name"))
    stats <- stats %>% dplyr::rename(bp = pos, gene = UCSC_RefGene_Name)
    gene <- stringr::str_split(stats$gene, pattern = ";", simplify = TRUE)
    get_unique_genes <- function(x) {
      unq <- unique(x)
      unq <- unq[unq != ""]
      paste(unq, sep = ";")
    }
    gene2 <- as.character(apply(gene, 1, FUN = get_unique_genes))
    gene2 <- ifelse(gene2 == "character(0)", "", gene2)
    stats$gene <- gene2
  }
  stats <- stats %>% dplyr::mutate(chr = as.numeric(stringr::str_replace(chr, "chr", "")))
  # Update gene names

  stats.tmp <-  stats %>%
    # Compute chromosome size
    dplyr::group_by(chr) %>%
    dplyr::summarise(chr_len=max(bp)) %>%

    # Calculate cumulative position of each chromosome
    dplyr::mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
    dplyr::select(-chr_len) %>%

    # Add this info to the initial dataset
    dplyr::left_join(stats, ., by=c("chr")) %>%

    # Add a cumulative position of each SNP
    dplyr::arrange(chr, bp) %>%
    dplyr::mutate( bpcum = bp + tot)

  if(!is.null(hlight)) {
    stats.tmp <- stats.tmp %>% dplyr::mutate(is_highlight = ifelse(Probe %in% hlight, TRUE, FALSE))
  }
  if(!is.null(annotate)) {
    stats.tmp <- stats.tmp %>% dplyr::mutate(is_annotate = ifelse(p < threshold, TRUE, FALSE))
  }

  # get chromosome center positions for x-axis
  axisdf <- stats.tmp %>% dplyr::group_by(chr) %>% dplyr::summarize(center=(max(bpcum) + min(bpcum) ) / 2 )

  g <- ggplot(stats.tmp, aes(x = bpcum, y = -log10(p))) +
    # Show all points
    geom_point(aes(color = as.factor(chr)), alpha= 0.8, size = point_size) +
    scale_color_manual(values = rep(col, 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis

    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome", y = expression(-log[10](italic(p)))) +

    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(threshold), linetype = "dashed") +

    # Custom the theme:
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      text = element_text(size = text_size),
      axis.text.x = element_text(hjust = 1, angle = 90, size = (text_size - 3)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  ## Options
  if(!is.null(sugg)) {
    g <- g + geom_hline(yintercept = -log10(sugg), linetype="dashed")
  }
  if(!is.null(hlight)) {
    g <- g + geom_point(data = stats.tmp %>% dplyr::filter(is_highlight), color = "red", size = point_size, shape = 23, fill = "red")
  }
  if(annotate) {
    g <- g + geom_label_repel(data = stats.tmp %>% dplyr::filter(is_annotate),
                              aes(label = as.factor(gene), alpha=0.7), size=5, force=1.3)
  }
  g

}


## Adapted from: https://gist.github.com/slowkow/9041570

gg_qqplot <- function(pval, show_lambda = TRUE, cut = FALSE, threshold = NULL, ci = NULL, text_size = 9, point_size = 1) {
  n  <- length(pval)
  chisq <- qchisq(1 - pval,1)
  inflation <- median(chisq)/qchisq(0.5,1)
  df <- tibble(
    observed = -log10(sort(pval)),
    expected = -log10(ppoints(n)),
    clower   = if(!is.null(ci)) -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)) else 0,
    cupper   = if(!is.null(ci)) -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)) else 0
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))

  if(cut) {
    if(is.null(threshold)) {
      threshold <- -log10(0.05 / length(pval))
    }
    df <- df %>%
      dplyr::mutate(observed = ifelse(observed > threshold, threshold, observed))
  }

  g <-  ggplot(df) +
    geom_point(aes(expected, observed), alpha= 0.8, size = point_size) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_classic() +
    theme(text = element_text(size = text_size))
  if(!is.null(ci)) {
    g <- g +
      geom_line(aes(expected, cupper), linetype = 2) +
      geom_line(aes(expected, clower), linetype = 2)
  }

  if(show_lambda) {

    label <- sprintf("paste(lambda, \" = \", %s)", round(inflation ,3))
    g <- g + annotate("text", y = max(df$observed), x = 0.7,
                      label = label, parse = TRUE, size = 0.31 * text_size)
  }
  g
}

compare_tstats <- function(stats1, stats2, label1, label2, highlight = NULL, highlight_sig = TRUE,
                           threshold = NULL, text_size = 9, point_size = 1, show_legend = TRUE,
                           annotation = NULL) {
  
  # Add logp
  if(!("t" %in% colnames(stats1))) {
    stats1 <- stats1 %>% dplyr::mutate(logp = -log10(p), t = b / se)
  }
  
  if(!("t" %in% colnames(stats2))) {
    stats2 <- stats2 %>% dplyr::mutate(logp = -log10(p), t = b / se)
  }
  
  if(is.null(threshold)) {
    threshold <- 0.05 / nrow(stats1)
  }
  
  combined <- stats1[,c("Probe", "p", "b", "t")] %>%
    dplyr::left_join(stats2[,c("Probe","p", "b", "t")], by = "Probe") %>%
    dplyr::mutate(Sig = dplyr::case_when(
      p.x < threshold & p.y < threshold ~ "Both Sig",
      p.x > threshold & p.y < threshold ~ sprintf("%s Sig", label2),
      p.x < threshold & p.y > threshold ~ sprintf("%s Sig", label1),
      TRUE ~ "Not Sig"
    )) %>% dplyr::mutate(Sig = factor(Sig, levels=c("Both Sig", sprintf("%s Sig", label2),
                                                    sprintf("%s Sig", label1), "Not Sig")))
  
  if(!is.null(highlight)) {
    combined <- combined %>% dplyr::mutate(highlight = ifelse(Probe %in% highlight, TRUE, FALSE))
  }
  
  max <- max(c(combined$t.x, combined$t.y), na.rm=T)
  min <- min(c(combined$t.x, combined$t.y), na.rm=T)
  
  color_values <- c("black", "#BC3C29", "#0072B5", "#E18727")
  color_values_names <- c("Not Sig", "Both Sig", sprintf("%s Sig", label2), sprintf("%s Sig", label1))
  names(color_values) <- color_values_names
  
  if(highlight_sig) {
    g <- ggplot(combined, aes(x = t.x,
                              y = t.y, color = Sig)) +
      geom_hline(yintercept = 0, linetype = 'dashed', colour = 'darkgrey') +
      geom_vline(xintercept = 0, linetype = 'dashed', colour = 'darkgrey') +
      geom_point(alpha= 0.8, size = point_size) +
      geom_point(data = combined %>% filter(Sig != "Not Sig"), alpha = 0.8, size = point_size) +
      geom_abline(slope = 1, linetype = "dashed") +
      scale_color_manual(values = color_values) +
      labs(x = sprintf("T-statistic %s", label1),
           y = sprintf("T-statistic %s", label2)) +
      theme_classic() +
      xlim(min - 5 , max + 5) +
      ylim(min - 5 , max + 5) +
      theme(text = element_text(size=text_size),
            legend.position = if(show_legend) "right" else "none") +
      coord_fixed()
  } else if (!is.null(annotation)) {
    combined <- combined %>%
      dplyr::left_join(annotation, by = "Probe")
    g <- ggplot(combined, aes(x = t.x,
                              y = t.y)) +
      geom_hline(yintercept = 0, linetype = 'dashed', colour = 'darkgrey') +
      geom_vline(xintercept = 0, linetype = 'dashed', colour = 'darkgrey') +
      geom_point(data = combined %>% filter(is.na(Annotation)), alpha= 0.8, size = point_size) +
      geom_point(data = combined %>% filter(!is.na(Annotation)),
                 aes(color = Annotation, shape = Annotation),
                 alpha = 1, size = point_size) +
      colorblindr::scale_color_OkabeIto() +
      geom_abline(slope = 1, linetype = "dashed") +
      labs(x = sprintf("T-statistic %s", label1),
           y = sprintf("T-statistic %s", label2)) +
      theme_classic() +
      xlim(min - 5 , max + 5) +
      ylim(min - 5 , max + 5) +
      theme(text = element_text(size=text_size),
            legend.position = if(show_legend) "right" else "none") +
      coord_fixed()
    
  } else {
    g <- ggplot(combined, aes(x = t.x,
                              y = t.y)) +
      geom_hline(yintercept = 0, linetype = 'dashed', colour = 'darkgrey') +
      geom_vline(xintercept = 0, linetype = 'dashed', colour = 'darkgrey') +
      geom_point(alpha= 0.8, size = point_size) +
      geom_abline(slope = 1, linetype = "dashed") +
      labs(x = sprintf("T-statistic %s", label1),
           y = sprintf("T-statistic %s", label2)) +
      theme_classic() +
      xlim(min - 5 , max + 5) +
      ylim(min - 5 , max + 5) +
      theme(text = element_text(size=text_size),
            legend.position = if(show_legend) "right" else "none") +
      coord_fixed()
  }
  
  
  if(!is.null(highlight)) {
    g <- g + geom_point(data = combined %>% dplyr::filter(highlight),
                        size = point_size, color = "red", shape = 8)
  }
  
  g
}


plot_pvalues_DMPs <- function(cpg, windowsize, DMPs,  correlation_plot = FALSE , nearby_probes = NULL, threshold = 1e-7,
                              significance_line=FALSE, show_legend = TRUE, text_size = 9, point_size = 1) {

  dmp <- DMPs[cpg,]
  # Find nearby CpGs
  cpgs <- subsetByOverlaps(DMPs, dmp, maxgap = round(windowsize/2), ignore.strand = TRUE)
  strand(cpgs) <- "*"
  cpgs <- sort(cpgs)

  # Create a tibble
  cpgs <- tibble(Probe = names(cpgs), Chr = as.character(seqnames(cpgs)),
                 p = cpgs$p,
                 logp = cpgs$logp, Pos = start(cpgs))

  # Add Significance
  cpgs <- cpgs %>% mutate(Sig = ifelse(p < threshold, TRUE, FALSE))
  max <- max(cpgs$logp)

  if(!significance_line) {
    plot <- ggplot(cpgs, aes(x = Pos, y = logp, color = Sig)) +
      geom_point(size = point_size) +
      xlab(cpgs$Chr[1]) +
      ylab(expression(-log[10](P))) +
      ggtitle(cpg) +
      scale_x_continuous(breaks = c(start(dmp) - round(windowsize/2) + 0.025 * windowsize,start(dmp), start(dmp) + round(windowsize/2) - 0.025 * windowsize),
                         labels = c(sprintf('- %s bp', round(windowsize/2)), as.character(start(dmp)), sprintf('+ %s bp', round(windowsize/2))),
                         limits = c(start(dmp) - round(windowsize/2), start(dmp) + round(windowsize/2) )) +
      scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
      ylim(0, max) +
      theme_classic()  +
      theme(text = element_text(size=text_size),
            axis.text = element_text(size=text_size),
            legend.position="none")

  } else {
    plot <- ggplot(cpgs, aes(x = Pos, y = logp, color = Sig)) +
      geom_point(size = point_size) +
      geom_hline(yintercept = -log10(threshold), linetype = 'dashed', color = 'red') +
      xlab(cpgs$Chr[1]) +
      ylab(expression(-log[10](P))) +
      ggtitle(sprintf('%s - %s', cpg, dmp$Gene)) +
      scale_x_continuous(breaks = c(start(dmp) - round(windowsize/2) + 0.025 * windowsize,start(dmp), start(dmp) + round(windowsize/2) - 0.025 * windowsize),
                         labels = c(sprintf('- %s bp', round(windowsize/2)), as.character(start(dmp)), sprintf('+ %s bp', round(windowsize/2))),
                         limits = c(start(dmp) - round(windowsize/2), start(dmp) + round(windowsize/2) )) +
      scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
      ylim(0, max) +
      theme_classic()  +
      theme(text = element_text(size=text_size),
            axis.text = element_text(size=text_size),
            legend.position="none")

  }

  ## Add correlation plot
  if(correlation_plot) {
    nearby.probes.tmp <- nearby_probes[[cpg]]
    beta <- nearby.probes.tmp$beta[cpgs$Probe,samplesheet_smallbatch$Sample_Name, drop=FALSE]
    cors <- cor(t(beta))
    cors <- melt(cors, varnames = c('Probe', 'Probe2'))

    labels <- cpgs$Probe
    labels <- ifelse(labels == cpg, cpg, "")
    heatmap <- cors %>% mutate(Probe = factor(Probe, levels = cpgs$Probe),
                               Probe2 = factor(Probe2, levels = unique(Probe))) %>%
      ggplot(aes(x = Probe, y = Probe2,  fill = value)) +
      geom_tile()  +
      scale_x_discrete(labels = labels) +
      scale_y_discrete(labels = labels) +
      # geom_text(aes(Probe, Probe2, label = round(value,2)), color = "black", size = 4) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0, limit = c(-1,1), space = "Lab",
                           name="Pearson\nCorrelation") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    return(list(plot, heatmap))
  } else {
    return(plot)
  }

}
