t2_plot <- function(p2,
                    show_ID_VAR = FALSE,
                    missing_legend_label = "Missing") {

  # Generate a Hotelling T-square plot

  # Args:
  #    p2: A list returned from function modP1
  #    show_ID_VAR: Defaulted to FALSE which only shows obs of UniqueID
  #    missing_legend_label: The label in legend for missing values

  # Returns:
  #    A Hotelling T-square plot

  new.scores <- p2$scores
  in.limit <- p2$UCL
  new.contributions <- p2$contributions
  new.OOC <- new.contributions[(OOC == "Y") %in% TRUE]
  new.NA <- new.contributions[is.na(OOC)]
  ID_VAR <- p2$id.var

  in.text <- paste("Number of OOCs =", nrow(new.OOC), sep = " ")
  if (nrow(new.NA) > 0) {
    in.text <- paste(in.text, ", Number of NAs = ",
      nrow(new.NA),
      sep = " "
    )
  }

  in.t2 <- new.scores[, c("T2", "OOC", ID_VAR), with = FALSE]

  # replace missing values with "M"
  in.t2[is.na(OOC), OOC := "M"]


  batch_size <- nrow(in.t2)

  if (show_ID_VAR == FALSE) {
    p.t2 <- ggplot(in.t2, aes(y = T2, x = seq(1, nrow(in.t2))), 
                   color = OOC) +
      theme_bw() +
      theme(
        axis.line = element_line(color = "gray"),
        panel.grid.major = element_blank()
      ) +
      xlab("obs") +
      theme(legend.position = "none")
  } else {
    p.t2 <- ggplot(in.t2, aes_string(x = ID_VAR, y = "T2", group = 1), 
                   color = OOC) +
      theme_bw() +
      theme(
        axis.line = element_line(color = "gray"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
      ) +
      xlab(ID_VAR)
  }

  p.t2 <- p.t2 +
    geom_point(aes(color = OOC)) +
    geom_hline(aes(yintercept = in.limit),
      color = "grey",
      linetype = "dashed", size = 1
    ) +
    ylab(expression(T^{
      2
    })) +
    ggtitle(expression(paste("Hotelling", T^2, "plot", sep = " "))) +
    annotate("text", x = 0.2 * batch_size, y = max(in.t2$T2) * 0.9, 
             label = in.text) +
    scale_colour_manual(
      name = "OOC",
      labels = c("Y" = "Y", "N" = "N", "M" = missing_legend_label),
      values = c("Y" = "red", "N" = "black", "M" = "orange")
    )


  return(p.t2)
}
