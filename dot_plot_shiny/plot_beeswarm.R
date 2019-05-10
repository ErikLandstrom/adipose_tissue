plot_beeswarm <- function(tb, gene_name) {
  tb %>% 
    dplyr::filter(name == gene_name) %>% 
    ggplot() +
    geom_beeswarm(aes(condition, intensity, color = missing),
                  size = 4,
                  groupOnX = TRUE) +
    labs(
      x = "Group",
      y = "log2(LFQ)",
      title = gene_name
    ) +
    scale_color_manual(name = "Imputed",
                       breaks = c("aObserved", "Missing"),
                       labels = c("Observed", "Missing"),
                       values = c("blue", "red")) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          plot.title = element_text(size = 18),
          axis.text = element_text(size = 12),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16))
}
