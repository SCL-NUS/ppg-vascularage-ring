# ============================================================
# Function: plot_per_fold_scatter
# Description:
#    Create scatter plot with per-fold correlation annotations

#created by Gizem @ 18 June
# ============================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)

# Define function
plot_per_fold_scatter <- function(data, xvar, yvar, fold_col,
                                  # x_label = xvar,
                                  # y_label = yvar,
                                  xlim = c(20, 70), ylim = c(20, 70)) {
  
  # Capture original variable names for axis labels
  x_label <- xvar
  y_label <- yvar

  # Dynamically extract the selected variables
  data <- data %>%
    mutate(
      xvar_val = .data[[xvar]],
      yvar_val = .data[[yvar]],
      Fold     = as.factor(.data[[fold_col]])   # <- THIS is the key line
      
    )
  
  # Calculate Pearson correlation coefficient per Fold
  r_squared_table <- data %>%
    group_by(Fold) %>%
    summarise(
      r = cor(xvar_val, yvar_val, method = "pearson"),
      .groups = 'drop'
    ) %>%
    arrange(Fold) %>%
    mutate(y_pos = ylim[2]+1 - (row_number() - 1) * 2)
  
  # Generate color palette for Folds
  folds <- sort(unique(data$Fold))
  n_folds <- length(folds)
  
  if (n_folds > 12) {
    stop("Paired palette can only handle up to 12 folds.")
  }
  
  my_palette <- brewer.pal(n_folds, "Paired")
  fold_colors <- data.frame(Fold = folds, color = my_palette)
  
  # Merge colors into R table for annotation coloring
  r_squared_table <- r_squared_table %>%
    left_join(fold_colors, by = "Fold")
  
  # Plot
  p <- ggplot(data, aes(x = xvar_val, y = yvar_val, color = Fold)) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey", linewidth = 0.5) +
    geom_point(size = 2, alpha = 0.6) + 
    scale_color_manual(values = setNames(my_palette, folds)) +
    labs(color = "Fold",x = x_label, y = y_label) +
    coord_cartesian(xlim = xlim, ylim = ylim) + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
    theme_minimal() +
    theme(
      # draw only x and y axis lines
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      # remove any remaining box/border
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      # remove grid
      panel.grid = element_blank(),
      # text formatting
      text = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    geom_text(
      data = r_squared_table,
      #aes(x = xlim[1] - 1, y = y_pos, label = paste0(Fold, ": R = ", sprintf("%.2f", r))),
      aes(x = xlim[1] - 1, y = y_pos, label = paste0("R = ", sprintf("%.2f", r))),
      
      color = r_squared_table$color,
      hjust = 0, vjust = 1, size = 4.5, inherit.aes = FALSE
    )
  
  return(p)
}



# 
# # Create mapping table between Fold and colors: to use in text annotations
# my_palette <- brewer.pal(length(unique(pred_CNN$Fold)), "Paired")
# fold_colors <- data.frame(
#   Fold = sort(unique(pred_CNN$Fold)),
#   color = my_palette
# )
# 
# # Calculate r and rsquared values for each fold:
# rm(r_squared_table)
# r_squared_table <- pred_CNN %>%
#   group_by(Fold) %>%
#   summarise(
#     r = cor(Chrono_Age, Estimated_Age, method = "pearson"),
#     r_squared = r^2
#   ) 
# # Add colors to your r_squared_table
# r_squared_table <- r_squared_table %>%
#   arrange(Fold) %>%
#   left_join(fold_colors, by = "Fold") %>%
#   mutate(y_pos = 71 - (row_number() - 1) * 2)
# 
# # Plot
# ggplot(pred_CNN, aes(x = Chrono_Age, y = Estimated_Age, color = Fold)) + 
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey", linewidth = 0.5) +
#   geom_point(size=2, alpha=0.6) + 
#   scale_color_brewer(palette = "Paired") +
#   labs(color = "Fold") +
#   coord_cartesian(xlim = c(20, 70), ylim = c(20, 70)) + 
#   geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .9) +
#   theme_minimal() +
#   theme(
#     axis.line = element_line(), 
#     text = element_text(size = 22),
#     axis.text.y=  element_text(size = 22),
#     axis.text.x=  element_text(size = 22),
#     plot.margin = unit(c(0, 0, 0, 0), "cm"),legend.position = "none",
#     panel.background = element_rect(fill = "white", colour = "grey50")
#   ) +
#   # Add R² labels for each Fold
#   geom_text(data = r_squared_table,
#             #aes(x = 22, y = y_pos, label = paste0(Fold, ": R² = ", sprintf("%.2f", r_squared))),
#             aes(x = 19, y = y_pos, label = paste0(Fold, ": R = ", sprintf("%.2f", r))),
#             color = r_squared_table$color,
#             hjust = 0, vjust = 1, size = 5)
