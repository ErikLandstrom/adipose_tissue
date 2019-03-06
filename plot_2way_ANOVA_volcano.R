### plot_2way_ANOVA_volcano
### Author: Erik Ländström
### Date: 190306


# Description -------------------------------------------------------------

# Plots a volcano plot for data generated with prepare_anova_without_int_for_volcano_plot.R

# Arguments ---------------------------------------------------------------

# tb = output from prepare_anova_without_int_for_volcano_plot

# Function ----------------------------------------------------------------

plot_2way_ANOVA_volcano <- function (tb) {
  tb %>% 
    ggplot(aes(lf2c, -log10(p_value), color = significant)) +
    geom_point() +
    scale_color_manual("More abundant in:",
                       values = c("black", "#00BFC4", "#F8666D", "#7CAE00", "#C77CFF"),
                       breaks = c(2, 3, 4, 5),
                       labels = c("MIDY", "WT", "Adipose, mesenterial", "Adipose, subcutaneous")) +
    facet_grid(~ group, scales = "free_x") +
    theme(panel.background = element_blank(),
          axis.line = element_line("black"))
}