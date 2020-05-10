contribution_plot <- function(p2,
                              contribution_rank = 0,
                              clustering_distance = "euclidean",
                              clustering_method = "ward.D",
                              color1sd = "black",
                              color2sd = "orange",
                              color3sd = "red",
                              vjust = -1,
                              hjust = -1) {

  # Generate two plots related to contributions

  # Args:
  #    p2: A list returned from function modP1
  #    contribution_rank: A flag for contribution plot. Defaulted to 0.
  #                       Use only OOCs for contributions rank plot when
  #                       it's set to 0. Use all for contributions rank plot
  #                       when it's set to 1. The variables on the X-axis of
  #                       heatmap are sorted by contributions with all rows when
  #                       it's set to 1. Plot both contributions rank plots
  #                       (for all and OOCs only) and the contribution heatmap
  #                       for OOCs when it's set to 2. The variables on the
  #                       X-axis of heatmap are sorted by contributions
  #                       with OOCs only when it's set to 2.
  #    clustering_distance: Distance metric for clustering observations
  #                        on the contribution heatmap
  #    clustering_method: Method for clustering observations
  #                      on the contribution heatmap
  #    color1sd: The color of line and text for 1sd on
  #              Rank of Variables Contributions plot
  #    color2sd: The color of line and text for 2sd on . Defaulted to black
  #              Rank of Variables Contributions plot. Defaulted to orange
  #    color3sd: The color of line and text for 3sd on
  #              Rank of Variables Contributions plot. Defaulted to red
  #    vjust: vjust of text for sd lines on Rank of Variables Contributions plot
  #    hjust jjust of text for sd lines on Rank of Variables Contributions plot

  # Returns:
  #    two contribution plots
  #    contribution.plot.data: A matrix used for the contribution heatmap
  #    contribution.plot.break: A break vector used for the contribution heatmap
  #    message: description for warning messages


  message <- "The program ran successfully."

  in.contribution <- p2$contributions[OOC == "Y"]

  new.scores <- p2$scores
  in.id_var <- p2$id.var
  in.vars <- p2$model.var
  all.mean <- p2$ic.mean
  all.sd <- p2$ic.sd

  if (contribution_rank %in% c(1, 2)) {

    # contribution rank plot for all
    all.smry <- p2$all_rank_smry
    all.pop_mean <- p2$all_rank_pop_mean
    # all.pop_sd <- p2$all_rank_pop_sd
    all.true_sd <- p2$all_rank_true_sd

    all.one_sigma <- all.pop_mean + all.true_sd
    all.two_sigma <- all.pop_mean + 2 * all.true_sd
    all.three_sigma <- all.pop_mean + 3 * all.true_sd

    if (!data.table::is.data.table(all.smry)) data.table::setDT(all.smry)

    all.smry[, c := ifelse(rank_group_mean > all.three_sigma, "3sd+",
      ifelse(rank_group_mean > all.two_sigma, "2sd+",
        ifelse(rank_group_mean > all.one_sigma, "1sd+", "0sd+")
      )
    )]

    all.smry <- all.smry[order(-rank_group_mean)]
    all.smry <- transform(all.smry, variable = reorder(variable, 
                                                       rank_group_mean))


    mypallet <- "Blues" # this pallet has mamximum 9 colors
    ordered.variables <- as.character(all.smry$variable)

    # plot0 is contribution rank plot for all
    plot0 <- ggplot(data = all.smry, aes(x = variable, y = rank_group_mean, 
                                         fill = factor(c))) +
      geom_bar(stat = "identity") +
      geom_hline(aes(yintercept = all.one_sigma), color = color1sd, 
                 linetype = "dashed", size = 1) +
      geom_text(aes(0, all.one_sigma, label = paste("1 Standard Deviation =",
        round(all.one_sigma, digits = 0),
        sep = " "
      )),
      size = 3, angle = 90, vjust = vjust, hjust = hjust, color = color1sd
      ) +
      geom_hline(aes(yintercept = all.two_sigma), color = color2sd, 
                 linetype = "dashed", size = 1) +
      geom_text(aes(0, all.two_sigma, label = paste("2 Standard Deviation =",
        round(all.two_sigma, digits = 0),
        sep = " "
      )),
      size = 3, angle = 90, vjust = vjust, hjust = hjust, color = color2sd
      ) +
      geom_hline(aes(yintercept = all.three_sigma), color = color3sd, 
                 linetype = "dashed", size = 1) +
      geom_text(aes(0, all.three_sigma, label = paste("3 Standard Deviation =",
        round(all.three_sigma, digits = 0),
        sep = " "
      )),
      size = 3, angle = 90, vjust = vjust, hjust = hjust, color = color3sd
      ) +
      coord_flip() +
      scale_fill_brewer(palette = mypallet) +
      ggtitle("Rank of Variables Contributions") +
      theme_bw() +
      theme(
        axis.line = element_line(color = "gray"),
        panel.grid.major = element_blank(),
        legend.position = "none",
        text = element_text(size = 10)
      )
  }

  if (nrow(in.contribution) > 0) {
    in.scores <- new.scores[OOC == "Y", c(in.id_var, "T2"), with = FALSE]
    in.contribution[in.scores, T2 := i.T2, on = in.id_var]
    raw.dat <- data.table::copy(in.contribution)
    n <- nrow(in.contribution)
  }

  if (contribution_rank %in% c(0, 2)) {

    # contribution rank plot for OOC
    if (nrow(in.contribution) > 0) {
      smry <- p2$OOC_rank_smry
      pop_mean <- p2$OOC_rank_pop_mean
      # pop_sd <- p2$OOC_rank_pop_sd
      true_sd <- p2$OOC_rank_true_sd


      one_sigma <- pop_mean + true_sd
      two_sigma <- pop_mean + 2 * true_sd
      three_sigma <- pop_mean + 3 * true_sd

      smry[, c := ifelse(rank_group_mean > three_sigma, "3sd+",
        ifelse(rank_group_mean > two_sigma, "2sd+",
          ifelse(rank_group_mean > one_sigma, "1sd+", "0sd+")
        )
      )]

      smry <- smry[order(-rank_group_mean)]
      smry <- transform(smry, variable = reorder(variable, rank_group_mean))


      mypallet <- "Blues" # this pallet has mamximum 9 colors
      ordered.variables <- as.character(smry$variable)

      # plot1 is contribution rank plot for OOC
      plot1 <- ggplot(data = smry, aes(x = variable, y = rank_group_mean, 
                                       fill = factor(c))) +
        geom_bar(stat = "identity") +
        geom_hline(aes(yintercept = one_sigma), color = color1sd, 
                   linetype = "dashed", size = 1) +
        geom_text(aes(0, one_sigma, label = paste("1 Standard Deviation =", 
                                                  round(one_sigma, digits = 0), 
                                                  sep = " ")),
          size = 3, angle = 90, vjust = vjust, hjust = hjust, color = color1sd
        ) +
        geom_hline(aes(yintercept = two_sigma), color = color2sd, 
                   linetype = "dashed", size = 1) +
        geom_text(aes(0, two_sigma, label = paste("2 Standard Deviation =", 
                                                  round(two_sigma, digits = 0), 
                                                  sep = " ")),
          size = 3, angle = 90, vjust = vjust, hjust = hjust, color = color2sd
        ) +
        geom_hline(aes(yintercept = three_sigma), color = color3sd, 
                   linetype = "dashed", size = 1) +
        geom_text(aes(0, three_sigma, label = paste("3 Standard Deviation =", 
                                                    round(three_sigma, 
                                                          digits = 0), 
                                                    sep = " ")),
          size = 3, angle = 90, vjust = vjust, hjust = hjust, color = color3sd
        ) +
        coord_flip() +
        scale_fill_brewer(palette = mypallet) +
        ggtitle("Rank of Variables Contributions") +
        theme_bw() +
        theme(
          axis.line = element_line(color = "gray"),
          panel.grid.major = element_blank(),
          legend.position = "none",
          text = element_text(size = 10)
        )
      # theme(legend.position="none")
    } else {
      message <- "No contribution rank plot is generatated due to 0 OOC"
      warning(message)
    }
  }



  if (nrow(in.contribution) > 0) {

    # plot 3: individual scaled contribution heatmap
    all.raw <- as.vector(as.matrix(raw.dat[, c(in.vars), with = FALSE]))

    if (is.na(all.mean) | is.na(all.sd)) {
      message <- "No contribution heatmap is generated due to 0 in control"
      warning(message)

      if (contribution_rank == 0) { # contribution rank plot for OOC only
        return(list(
          plot1 = plot1, message = message,
          contribution.plot.data = NULL, contribution.plot.break = NULL
        ))
      } else if (contribution_rank == 1) { # contribution rank plot for ALL only
        return(list(
          plot0 = plot0, message = message,
          contribution.plot.data = NULL, contribution.plot.break = NULL
        ))
      } else if (contribution_rank == 2) { # contribution rank plot for both
        return(list(
          plot0 = plot0, plot1 = plot1, message = message,
          contribution.plot.data = NULL, contribution.plot.break = NULL
        ))
      }
    }

    ymax <- (max(all.raw) - all.mean) / all.sd
    ymin <- (min(all.raw) - all.mean) / all.sd

    # b <- c(ymin,1,2,3,ymax)
    ### improve color scales in the heatmap
    if (ymax > 3) {
      b <- c(ymin, 
             3, 
             (ymax - 3) / 8, 
             (ymax - 3) / 4, 
             (ymax - 3) / 8 * 3, 
             (ymax - 3) / 2, 
             (ymax - 3) / 8 * 5, 
             (ymax - 3) / 4 * 3, 
             (ymax - 3) / 8 * 7, 
             ymax)
    } else {
      b <- c(ymin, 1, 2, 3)
    }
    ####################

    graph.dat.3 <- data.table::copy(raw.dat)
    graph.dat.3[, (in.vars) := lapply(.SD, function(x) (x - all.mean) / all.sd), 
                .SDcols = in.vars]
    graph.dat.3 <- graph.dat.3[order(-T2, graph.dat.3[[in.id_var]])]

    ordered.sns <- graph.dat.3[[in.id_var]] # return a vector

    colsToDelete <- c("T2", "OOC", in.id_var)
    graph.dat.3[, (colsToDelete) := NULL]

    data.table::setcolorder(graph.dat.3, ordered.variables)


    mat.data.3 <- data.matrix(graph.dat.3)
    rownames(mat.data.3) <- ordered.sns
    colnames(mat.data.3) <- ordered.variables

    p <- length(in.vars)
    width <- max(4, ceiling(p / 3))
    height <- max(4, ceiling(n / 6))

    # plot2 is contribution heatmap for OOC

    if (nrow(mat.data.3) > 2) {
      plot2 <- pheatmap(mat.data.3,
        col = brewer.pal(length(b) - 1, mypallet),
        main = "Globally Standardized Contribution",
        height = height, width = width,
        breaks = b, scale = "none", cluster_rows = T, cluster_cols = F,
        clustering_distance_rows = clustering_distance, 
        clustering_distance_cols = clustering_distance,
        clustering_method = clustering_method
      )
    } else {
      plot2 <- pheatmap(mat.data.3,
        col = brewer.pal(length(b) - 1, mypallet),
        main = "Globally Standardized Contribution",
        height = height, width = width,
        breaks = b, scale = "none", cluster_rows = F, cluster_cols = F
      )
    }

    if (contribution_rank == 0) { # contribution rank plot for OOC only
      return(list(
        plot1 = plot1, plot2 = plot2, message = message,
        contribution.plot.data = mat.data.3, contribution.plot.break = b
      ))
    } else if (contribution_rank == 1) { # contribution rank plot for ALL only
      return(list(
        plot0 = plot0, plot2 = plot2, message = message,
        contribution.plot.data = mat.data.3, contribution.plot.break = b
      ))
    } else if (contribution_rank == 2) { # contribution rank plot for both
      return(list(
        plot0 = plot0, plot1 = plot1, plot2 = plot2, message = message,
        contribution.plot.data = mat.data.3, contribution.plot.break = b
      ))
    }
  } else {
    message <- "No contribution rank plotis generatated due to 0 OOC"
    warning(message)

    if (contribution_rank == 0) { # contribution rank plot for OOC only
      return(list(
        message = message,
        contribution.plot.data = NULL, contribution.plot.break = NULL
      ))
    } else if (contribution_rank == 1) { # contribution rank plot for ALL only
      return(list(
        plot0 = plot0, message = message,
        contribution.plot.data = NULL, contribution.plot.break = NULL
      ))
    } else if (contribution_rank == 2) { # contribution rank plot for both
      return(list(
        plot0 = plot0, message = message,
        contribution.plot.data = NULL, contribution.plot.break = NULL
      ))
    }
  }
}
