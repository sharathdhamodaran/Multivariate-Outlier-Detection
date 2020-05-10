mod_p2 <- function(p1,
                   in.score_data,
                   na.rm = TRUE,
                   na.T2 = NA,
                   batch_alpha = 0.01,
                   batch_beta = 0.01) {

  # Phase II - Score on new data and get contribution

  # Args:
  #    p1: A list returned from function modP1
  #    in.score_data: A data frame or data table for scoring
  #    na.rm: Whether NA values should be stripped before computation proceeds
  #    na.T2: A user specified T2 value for rows with missing parameters
  #           It's used only when na.rm is set to FALSE.  If left as NA
  #           while na.rm is set to FALSE, the Hotelling T-squres are calculated
  #           for all individuals/rows, ingoring the missing columns for
  #           individuals/rows with missing parametric values.
  #    batch_alpha: type I error for batch OOC decision rule
  #    batch_beta: type II error for batch OOC decision rule

  # Returns:
  #    scores: A list with a score data table
  #    contributions: A contribution data table
  #    ic.mean: Mean of in-control objects' contributions
  #    ic.sd: Standar Deviation of in-control objects' contributions
  #    model.var: All variables used in modP1
  #    id.var: An id variable
  #    UCL: Upper control limit
  #    UCL_pct: Percentage used for UCL
  #    batch_OOC: A flag indicating whether the batch is OOC or not
  
  #    all_rank_smry: A datatable with all the parametric variables and their
  #                   mean contribution to OOC
  #    all_rank_pop_mean: the mean of all variables 
  #    all_rank_pop_sd: standard deviation of all variables
  #    all_rank_true_sd: expected standard deviation of all rank means of all
  #                      variables
  #    OOC_rank_smry: A datatable with all the parametric variables and their
  #                   mean contribution to OOC if there's any OOC
  #    OOC_rank_pop_mean: the mean of all variables of OOC individuals
  #    OOC_rank_pop_sd: standard deviation of all variables of OOC individuals
  #    OOC_rank_true_sd: expected standard deviation of all rank means of all
  #                      variables of OOC individuals
  #    big.contributors: a vector of characters containing the biggeste
  #                      contributors
  #    message = description of messages for warnings and execution
  
  message <- "The program ran successfully."

  in.S <- p1$S
  in.xbar <- p1$xbar
  in.vars <- p1$model.var

  in.UCL <- p1$UCL
  in.id_var <- p1$id.var

  weight <- p1$weight

  if (!data.table::is.data.table(in.score_data)) {
    in.score_data <- data.table::setDT(in.score_data)
  }

  if ((na.rm == FALSE) & is.na(na.T2)) {
    xdata_new <- in.score_data
  } else {
    xdata_new <- in.score_data[complete.cases(
      in.score_data[, in.vars, with = FALSE]
    )]
  }

  # store original order
  orig.order <- in.score_data[, list(row = list(.I)), by = in.id_var]

  if ((na.rm == FALSE) & (!is.na(na.T2))) {
    null_obs <- in.score_data[!complete.cases(
      in.score_data[, in.vars, with = FALSE]
    ),
    setdiff(colnames(in.score_data), in.vars),
    with = FALSE
    ]
  }

  xdata <- xdata_new[, in.vars, with = FALSE]

  xbar <- as.vector(in.xbar)
  deltaX_new <- as.matrix(xdata - xbar[col(xdata)])

  # make a copy to get index of NA values
  deltaX_orig <- deltaX_new

  # set NA to 0
  deltaX_new[is.na(deltaX_orig)] <- 0

  if (!is.null(p1$inv_S)) {
    inv_S <- p1$inv_S
  } else {
    inv_S <- tryCatch(solve(in.S), error = function(e) {
      return(NULL)
    })

    if (is.null(inv_S)) {
      inv_S <- tryCatch(MASS::ginv(in.S), error = function(e) {
        return(NULL)
      })
    }

    if (is.null(inv_S)) {
      message <- "Singular matrix. Remove some columns"
      warning(message)
      return(message = message)
    }
  }

  # calculate Hotelling T-square

  if (is.null(weight)) {
    T2.new <- colSums(t(deltaX_new %*% inv_S) * t(deltaX_new))
  } else {
    weighted.inv_S <- diag(weight) %*% inv_S

    T2.new <- colSums(t(deltaX_new %*% weighted.inv_S) * t(deltaX_new))
  }

  invisible(gc())

  scores <- cbind(xdata_new, T2 = T2.new)

  scores$OOC <- ifelse(scores$T2 > in.UCL, "Y", "N")

  # contributions
  p <- length(in.vars)

  for (pp in 1:p) {
    xdat_i <- data.table::copy(xdata)
    xdat_i <- xdat_i[, (pp) := NULL]

    xbar_i <- t(matrix(xbar[-pp]))
    deltaX_i <- as.matrix(xdat_i - xbar_i[col(xdat_i)])

    deltaX_i_orig <- deltaX_i

    # set NA to 0
    deltaX_i[is.na(deltaX_i_orig)] <- 0

    if (!is.null(p1$inv_S_i)) {
      inv_S_i <- p1$inv_S_i[[pp]]
    } else {
      S_i <- in.S[-pp, -pp]

      inv_S_i <- tryCatch(solve(S_i), error = function(e) {
        return(NULL)
      })

      if (is.null(inv_S_i)) {
        inv_S_i <- tryCatch(MASS::ginv(S_i), error = function(e) {
          return(NULL)
        })
      }

      if (is.null(inv_S_i)) {
        message <- "Singular matrix. Remove some columns"
        warning(message)
        return(message = message)
      }
    }

    if (is.null(weight)) {
      D2_i <- T2.new - colSums(t(deltaX_i %*% inv_S_i) * t(deltaX_i))
    } else {
      m_weight <- diag(weight)

      weight_i <- m_weight[-pp, -pp]

      weighted.inv_S_i <- weight_i %*% inv_S_i

      D2_i <- T2.new - colSums(t(deltaX_i %*% weighted.inv_S_i) * t(deltaX_i))
    }


    invisible(gc())

    if (pp == 1) {
      Contributions <- data.table::as.data.table(D2_i)
      data.table::setnames(Contributions, "D2_i", in.vars[pp])
    } else {
      Contributions[, (in.vars[pp]) := D2_i]
    }
  }

  # set the contributions of cells with missing value to 0
  Contributions[is.na(deltaX_orig)] <- 0

  Contributions[, (in.id_var) := xdata_new[, in.id_var, with = FALSE]]
  Contributions[, "OOC" := scores$OOC]

  new.raw <- as.vector(as.matrix(Contributions[(OOC == "N") %in% TRUE, 
                                               c(in.vars), with = FALSE]))
  new.mean <- mean(new.raw, na.rm = TRUE)
  new.sd <- sd(new.raw, na.rm = TRUE)

  new.OOC <- Contributions[(OOC == "Y") %in% TRUE]

  if (nrow(xdata_new) > 1) {
    sample_size_result <- mspc::unknown_dec_rule(
      batch_size = nrow(xdata_new),
      t2_threshold_quantile = 1 - p1$UCL_pct / 100,
      t2_outlier_quantile = p1$t2_outlier_quantile,
      alpha = batch_alpha,
      beta = batch_beta
    )

    if (nrow(new.OOC) < sample_size_result[[1]]) {
      batch_OOC <- "N"
    } else {
      batch_OOC <- "Y"
    }
  } else if (nrow(xdata_new) == 1) {
    if (nrow(new.OOC) == 0) {
      batch_OOC <- "N"
    } else {
      batch_OOC <- "Y"
    }
  }

  # add contributions for all
  all.scores <- scores[, c(in.id_var, "T2"), with = FALSE]
  all.contributions <- data.table::copy(Contributions)

  p <- length(in.vars)

  all.contributions[all.scores, T2 := i.T2, on = in.id_var]
  all.contributions[, sumx := rowSums(.SD, na.rm = T), .SDcols = in.vars]
  all.contributions[, (in.vars) := lapply(.SD, function(x) x / sumx), 
                    .SDcols = in.vars]
  # colsToDelete <- c("sumx","OOC")
  all.contributions[, sumx := NULL]

  n <- nrow(all.contributions)
  N <- n * p

  new.all.long <- reshape2::melt(all.contributions,
    id.vars = c(in.id_var, "T2", "OOC"), measures.vars = in.vars,
    variable.name = "variable", value.name = "contribution"
  )

  # calculate the actual group means (rank)
  all.pop_mean <- (N + 1) / 2 # same as group mean
  all.pop_sd <- sqrt((N^2 - 1) / (12 * n))
  all.true_sd <- sqrt((N + 1) * (p - 1) / 12) # sd of group mean

  data.table::setDT(new.all.long)
  new.all.long[, rank := rank(contribution)]
  all.smry <- new.all.long[, .(rank_group_mean = mean(rank)), 
                           by = variable][order(-rank_group_mean)]

  # add code here to get big contributors
  if (nrow(new.OOC) > 0) {

    n <- nrow(new.OOC)
    N <- n * p

    new.OOC.long <- new.all.long[OOC == "Y"]

    # calculate the actual group means (rank)
    pop_mean <- (N + 1) / 2 # same as group mean
    pop_sd <- sqrt((N^2 - 1) / (12 * n))
    true_sd <- sqrt((N + 1) * (p - 1) / 12) # sd of group mean

    # data.table::setDT(new.OOC.long)
    new.OOC.long[, rank := rank(contribution)]
    smry <- new.OOC.long[, .(rank_group_mean = mean(rank)), 
                         by = variable][order(-rank_group_mean)]
    big.contributors <- as.character(smry[rank_group_mean > pop_mean + 
                                            3 * true_sd]$variable)
  } else {
    smry <- NULL
    pop_mean <- NA
    pop_sd <- NA
    true_sd <- NA
    big.contributors <- NULL
  }

  ### add steps for handling NAs
  if ((na.rm == FALSE) & (!is.na(na.T2))) {
    if (nrow(null_obs) > 0) {
      null_obs[, c(in.vars, "OOC") := NA]
      null_obs[, T2 := na.T2]
      # scores <- rbind(scores,null_obs, stringsAsFactors = FALSE)
      scores <- data.table::rbindlist(list(scores, null_obs), use.names = TRUE)

      # merge to get original order
      scores <- merge(orig.order, scores, by = in.id_var, sort = FALSE)
      scores[, row := NULL]

      null_obs <- null_obs[, c(in.vars, in.id_var, "OOC"), with = FALSE]
      # Contributions <- rbind(Contributions,null_obs,stringsAsFactors = FALSE)
      Contributions <- data.table::rbindlist(list(Contributions, null_obs), 
                                             use.names = TRUE)

      # merge to get original order
      Contributions <- merge(orig.order, Contributions, by = in.id_var, 
                             sort = FALSE)
      Contributions[, row := NULL]
    }
  }

  return(list(
    scores = scores, contributions = Contributions,
    ic.mean = new.mean, ic.sd = new.sd,
    model.var = in.vars, id.var = in.id_var,
    UCL = in.UCL, UCL_pct = p1$UCL_pct,
    batch_OOC = batch_OOC,
    all_rank_smry = all.smry,
    all_rank_pop_mean = all.pop_mean,
    all_rank_pop_sd = all.pop_sd,
    all_rank_true_sd = all.true_sd,
    OOC_rank_smry = smry,
    OOC_rank_pop_mean = pop_mean,
    OOC_rank_pop_sd = pop_sd,
    OOC_rank_true_sd = true_sd,
    big.contributors = big.contributors,
    message = message
  ))
}
