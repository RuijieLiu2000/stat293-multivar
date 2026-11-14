# FIGURES AND TABLES

# Input: Metric summaries and destination folder.
# Output: Four PNG figures saved in the folder.
plot_results <- function(summary, out_dir) {
  factors <- c("innovation_scale", "dimension", "density",
               "contemporaneous_correlation", "serial_correlation", "measurement_noise")
  columns <- c("sigma", "d", "density", "cross_rho", "serial_rho", "measurement_sd")
  labels <- c("Innovation SD", "Number of variables d", "Nonzero fraction",
              "Contemporaneous correlation", "Serial error correlation", "Measurement error SD")
  methods <- c("multiVAR_FISTA", "adaptive_multiVAR", "individual_lasso", "pooled_lasso")
  method_labels <- c("Multi-VAR (FISTA)", "Adaptive multi-VAR", "Individual LASSO", "Pooled LASSO")
  colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  for (metric in c("relative_frobenius", "rmsfe_h1", "sensitivity", "specificity")) {
    png(file.path(out_dir, paste0(metric, ".png")), width = 1800, height = 1250, res = 150)
    layout(matrix(c(1:6, 7, 7, 7), nrow = 3, byrow = TRUE), heights = c(1, 1, .16))
    par(mar = c(4.1, 4.3, 2.5, 1), oma = c(2, 0, 1, 0), las = 1)
    for (i in seq_along(factors)) {
      z <- summary[summary$metric == metric & summary$factor %in% c("baseline", factors[i]), ]
      yr <- range(c(z$lower, z$upper, z$mean), na.rm = TRUE)
      y_label <- switch(metric, relative_frobenius = "Relative Frobenius error",
                        rmsfe_h1 = "One-step RMSFE", sensitivity = "Sensitivity", specificity = "Specificity")
      plot(range(z[[columns[i]]]), yr, type = "n", xlab = labels[i],
           ylab = y_label, main = labels[i], cex.main = 1)
      grid(col = "#EEEEEE")
      for (m in seq_along(methods)) {
        zz <- z[z$method == methods[m], ]
        zz <- zz[order(zz[[columns[i]]]), ]
        x <- zz[[columns[i]]]
        lines(x, zz$mean, col = colors[m], lwd = 1.8, type = "b", pch = 14 + m)
        ok <- is.finite(zz$lower) & zz$upper > zz$lower
        if (any(ok)) arrows(x[ok], zz$lower[ok], x[ok], zz$upper[ok],
                            angle = 90, code = 3, length = .035, col = colors[m])
      }
    }
    par(mar = rep(0, 4))
    plot.new()
    legend("center", legend = method_labels, col = colors, pch = 15:18,
                       lty = 1, horiz = TRUE, bty = "n", cex = .8)
    mtext("Points: replicate means. Bars: 95% Monte Carlo t intervals (pointwise).", outer = TRUE,
          side = 1, line = .5, cex = .8)
    dev.off()
  }
}

# Input: None.
# Output: Display labels indexed by method name.
method_labels <- function() {
  c(multiVAR_FISTA = "Multi-VAR", adaptive_multiVAR = "Adaptive",
    individual_lasso = "Individual", pooled_lasso = "Pooled")
}

# Input: Scenario table with factor names and values.
# Output: One display label per scenario.
scenario_labels <- function(scenarios) {
  labels <- c(baseline = "Baseline", innovation_scale = "Innovation SD",
              dimension = "Variables", density = "Density",
              contemporaneous_correlation = "Cross correlation",
              serial_correlation = "Serial correlation",
              measurement_noise = "Measurement SD", series_length = "Time points")
  ifelse(scenarios$factor == "baseline", "Baseline",
         paste(unname(labels[scenarios$factor]), scenarios$value))
}


# Input: Table, title, note, and destination path.
# Output: Writes the Markdown table and returns its path invisibly.
write_display_table <- function(table, title, note, path) {
  lines <- c(paste0("**", title, "**"), "", note, "",
             knitr::kable(table, format = "pipe", row.names = FALSE), "")
  writeLines(lines, path)
  invisible(path)
}


# Input: Metric summaries or contrasts, scenarios, metric, value column, and methods.
# Output: Display table of estimates and SEs by scenario and method.
metric_table <- function(values, scenarios, metric, value_column = "mean",
                         methods = names(method_labels())) {
  table <- data.frame(Condition = scenario_labels(scenarios))
  for (method in methods) {
    z <- values[values$metric == metric & values$method == method, ]
    z <- z[match(scenarios$scenario, z$scenario), ]
    table[[method_labels()[[method]]]] <- sprintf(
      "%.3f (%.3f)", z[[value_column]], z$mcse
    )
  }
  table
}


# Input: Replicate results, metric summaries, scenario settings, and destination folder.
# Output: Four figures and seven tables saved under the destination folder.
write_key_results <- function(raw, summary, scenarios, output_dir = "results") {
  figure_dir <- file.path(output_dir, "figures")
  table_dir <- file.path(output_dir, "tables")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  plot_results(summary, figure_dir)

  design <- data.frame(
    Condition = scenario_labels(scenarios), d = scenarios$d, T = scenarios$T,
    Density = scenarios$density, SD = scenarios$sigma,
    Cross = scenarios$cross_rho, Serial = scenarios$serial_rho,
    Measurement = scenarios$measurement_sd
  )
  write_display_table(
    design, "Simulation design",
    "All conditions: 10 subjects, common fraction 0.5, maximum spectral radius 0.8.",
    file.path(table_dir, "00_scenario_design.txt")
  )

  titles <- c(relative_frobenius = "Relative Frobenius error",
              rmsfe_h1 = "One-step forecast RMSE", sensitivity = "Sensitivity",
              specificity = "Specificity")
  note <- paste0(
    "Values are mean (Monte Carlo SE); ", length(unique(raw$replication)),
    " replications per condition. Adaptive = adaptive multi-VAR; ",
    "Individual/Pooled = individual/pooled LASSO."
  )
  for (i in seq_along(titles)) {
    metric <- names(titles)[i]
    write_display_table(
      metric_table(summary, scenarios, metric), titles[[metric]], note,
      file.path(table_dir, sprintf("%02d_%s.txt", i, metric))
    )
  }

  focus <- c("relative_frobenius", "rmsfe_h1")
  contrasts <- paired_contrasts(raw, focus)
  for (i in seq_along(focus)) {
    metric <- focus[i]
    write_display_table(
      metric_table(contrasts, scenarios, metric, "difference",
                   c("multiVAR_FISTA", "adaptive_multiVAR", "pooled_lasso")),
      paste("Paired difference:", titles[[metric]]),
      paste("Method minus individual LASSO, paired within replicate.",
            "Values are mean difference (Monte Carlo SE); negative favors the method."),
      file.path(table_dir, sprintf("%02d_paired_%s.txt", i + 4L, metric))
    )
  }
  invisible(NULL)
}
