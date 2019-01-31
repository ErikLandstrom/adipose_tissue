
plot_fdr_volcano <- function(tb, plot_title = "Volcano plot", break_vector, color_vector) {
  tb %>% 
    ggplot(aes(log2_difference, -log10(p_value), color = factor(significant))) +
    geom_point() +
    xlab(label = "Log2 fold change)") +
    ylab(label = "-Log10(p-value") +
    labs(
      title = plot_title
    ) +
    scale_color_manual(name = "More abundant in:",
                       breaks = break_vector,
                       values = color_vector) +
    theme(panel.background = element_blank(),
          axis.line = element_line("black"),
          legend.key = element_blank())
}
