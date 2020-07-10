manhattanPlot <- function(assocStats,
                          trait,
                          threshold,
                          colors = c("#91baff", "#3e619b"),
                          showSigMarkers = FALSE) {
  
  ## Get chromosome count
  chrs <- length(levels(assocStats$Chr))
  
  ## Filter data
  filtglmStats <- assocStats %>%
    dplyr::filter(.data$Trait == trait) %>%
    dplyr::mutate(highlight_flag = ifelse(-log10(p) >= threshold, T, F))
  
  ## Color data
  col <- c(
    "col1" = "#333333",
    "col2"   = "#D34747"
  )
  
  ## Plot components
  p <- ggplot2::ggplot(data = filtglmStats) +
    ggplot2::aes(x = .data$Pos, y = -log10(.data$p)) +
    ggplot2::geom_point(size = 0.8, na.rm = TRUE) +
    ggplot2::aes(color = .data$Chr) +
    ggplot2::geom_point(
      data = dplyr::filter(
        filtglmStats, .data$highlight_flag == T
      ),
      color = col[["col1"]],
      size = 2
    ) +
    ggplot2::scale_color_manual(
      values = rep(colors, chrs)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed") +
    ggplot2::xlab("Position") +
    ggplot2::ylab(bquote(~-log[10]~ '('*italic(p)*'-value)')) +
    ggplot2::geom_rug(sides = "b", colour = "black", alpha = 0.1) +
    ggplot2::ggtitle(label = paste("Trait:", trait)) +
    ggplot2::facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
    if (!showSigMarkers) {
      NULL
    } else {
      ggplot2::geom_text(
        ggplot2::aes(
          label = ifelse(
            test = -log10(p) >= threshold,
            yes = as.character(.data$Marker),
            no = ""
          )
        ),
        color = col[["col2"]],
        size = 3,
        hjust = "inward",
        vjust = "inward",
        na.rm = TRUE
      )
    }
  return(p)
}
