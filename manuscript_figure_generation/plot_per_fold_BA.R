
plot_per_fold_BA <- function(
    data,
    xvar,
    yvar,
    fold_col,
    xlim = c(20, 70),
    ylim = c(-30, 30),
    use_mean_x = FALSE,              # FALSE: delta ~ chrono (xvar); TRUE: delta ~ mean(x,y)
    delta_adjusted = FALSE,              # defauls FALSE: 
    annotate_x = NULL,               # where to place fold slope text; default uses xlim[1]
    annotate_y_top = NULL,           # top y for text; default uses ylim[2]
    annotate_line_step = 2           # vertical spacing between text lines
) {
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  
  # Axis labels (defaults to variable names)
  x_label <- xvar
  y_label <- yvar
  
  # Prepare data
  data2 <- data %>%
    mutate(
      xvar_val = .data[[xvar]],
      yvar_val = .data[[yvar]],
      Fold     = as.factor(.data[[fold_col]]),
      deltaAge = yvar_val - xvar_val,
      mean_x   = (xvar_val + yvar_val) / 2
    )
  
  # Choose x for BA regression
  data2 <- data2 %>%
    mutate(x_ba = if (use_mean_x) mean_x else xvar_val)
  
  # Bland–Altman summary (overall, as in your current code)
  ba_stats <- data2 %>%
    summarise(
      mean_bias = mean(deltaAge, na.rm = TRUE),
      sd_bias   = sd(deltaAge, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      upper_limit = mean_bias + 1.96 * sd_bias,
      lower_limit = mean_bias - 1.96 * sd_bias
    )
  
  mean_bias   <- ba_stats$mean_bias
  upper_limit <- ba_stats$upper_limit
  lower_limit <- ba_stats$lower_limit
  
  # Fold palette
  folds <- sort(unique(data2$Fold))
  n_folds <- length(folds)
  if (n_folds > 12) stop("Paired palette can only handle up to 12 folds.")
  my_palette <- brewer.pal(n_folds, "Paired")
  names(my_palette) <- folds
  
  # ---- Per-fold OLS slopes + 95% CI ----
  # slope CI via standard lm CI
  slope_table <- data2 %>%
    group_by(Fold) %>%
    do({
      df <- .
      fit <- lm(deltaAge ~ x_ba, data = df)
      
      ci <- suppressMessages(confint(fit, level = 0.95))
      # rownames: (Intercept), x_ba
      tibble(
        slope      = coef(fit)[["x_ba"]],
        slope_lwr  = ci["x_ba", 1],
        slope_upr  = ci["x_ba", 2]
      )
    }) %>%
    ungroup() %>%
    arrange(Fold) %>%
    mutate(
      # text placement
      y_pos = (if (is.null(annotate_y_top)) ylim[2] else annotate_y_top) -
        (row_number() - 1) * annotate_line_step,
      x_pos = (if (is.null(annotate_x)) xlim[1] else annotate_x),
      # label = paste0(
      #   Fold, ": slope = ", sprintf("%.2f", round(slope,2)),
      #   " [", sprintf("%.2f", round(slope_lwr,2)), ", ", sprintf("%.2f", round(slope_upr,2)), "]"
      # ),
      label = paste0(
        "slope [CI] = ", sprintf("%.2f", round(slope,2)),
        " [", sprintf("%.2f", round(slope_lwr,2)), ", ", sprintf("%.2f", round(slope_upr,2)), "]"
      ),
      color = my_palette[as.character(Fold)]
    )
  
  # ---- Plot ----
  offset <- 2
  
  p <- ggplot(data2, aes(x = x_ba, y = deltaAge, color = Fold)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = my_palette) +
    
    # OLS line per fold
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    
    # Overall BA lines (mean bias & LoA)
    geom_hline(yintercept = mean_bias,  color = "azure4", linewidth = 0.4) +
    geom_hline(yintercept = upper_limit, linetype = "dashed", color = "azure4", linewidth = 0.4) +
    geom_hline(yintercept = lower_limit, linetype = "dashed", color = "azure4", linewidth = 0.4) +
    
    # numeric annotations for BA lines
    annotate("text", x = xlim[2] - 4, y = mean_bias + offset,
             label = paste0(round(mean_bias, 2)), hjust = 0, size = 5) +
    annotate("text", x = xlim[2] - 4, y = upper_limit + offset,
             label = paste0(round(upper_limit, 2)), hjust = 0, size = 5) +
    annotate("text", x = xlim[2] - 4, y = lower_limit + offset,
             label = paste0(round(lower_limit, 2)), hjust = 0, size = 5) +
    
    # Fold slope + CI "legend" at the top (colored)
    geom_text(
      data = slope_table,
      aes(x = x_pos, y = y_pos, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1, size = 4.5,
      color = slope_table$color
    ) +
    
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(
      x = if (use_mean_x) "Mean age ((Chrono + Estimated)/2)" else x_label,
      y = if (delta_adjusted) "ΔAge_adjusted "
        else "ΔAge (Estimated − Chronological)",
      color = "Fold"
    ) +
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
    )
  
  return(p)
}
