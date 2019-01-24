### plot_volcano_pvalue
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Plots an volcano from data generated with multiple_ttests.R.


# Arguments ---------------------------------------------------------------

# tb = output from multiple_ttests.R
# label_1 = character, more abundant in " "
# label_2 = character, more abundant in " "
# plot_title = plot title


# Function ----------------------------------------------------------------

plot_volcano_pvalue <- function(tb, label_1, label_2,
                                plot_title = "Volcano plot") {
  
  # Volcano_plot
  tb %>% 
    mutate(
      enriched_1 = as.numeric(log2_difference > 0.585 & p_value < 0.05),
      enriched_2 = as.numeric(log2_difference < -0.585 & p_value < 0.05),
      significant = as.factor((enriched_1 + enriched_2 * 2) + 1)
    ) %>%
    ggplot(aes(log2_difference, -log10(p_value), color = significant)) +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = log2(1.5)) +
    geom_vline(xintercept = -log2(1.5)) +
    geom_point() +
    scale_color_manual("More abundant in:",
                       values = c("black", "#00BFC4", "#F8666D"),
                       breaks = c(2, 3),
                       labels = c(label_1, label_2)
    ) +
    labs(
      title = plot_title
    ) +
    xlab(bquote(~ log[2] ~ " fold change")) +
    ylab(bquote("-" ~ log[10] ~ "(p-value)")) +
    #  title = plot_title
    # ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line("black")
    )
}

